use anyhow::{Context, Result, anyhow};
use std::cmp::{max, min};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// A genomic interval (1-based, inclusive on both ends).
#[derive(Debug, Clone)]
pub struct Region {
    pub name: Option<String>,
    pub chr: String,
    pub start: u64,
    pub end: u64,
}

/// Parse comma-separated inline region strings.
pub fn parse_regions_inline(s: &str, flank: Option<u64>) -> Result<Vec<Region>> {
    s.split(',')
        .filter(|t| !t.trim().is_empty())
        .map(|t| parse_region_str(t.trim(), flank))
        .collect()
}

/// Parse a file with one region string per line.
/// If `path` is "-", reads from stdin instead.
pub fn parse_regions_list(path: &Path, flank: Option<u64>) -> Result<Vec<Region>> {
    let reader: Box<dyn BufRead> = if path == Path::new("-") {
        Box::new(BufReader::new(std::io::stdin()))
    } else {
        let f = File::open(path)
            .with_context(|| format!("Cannot open list file: {}", path.display()))?;
        Box::new(BufReader::new(f))
    };
    let source = path.display().to_string();
    reader
        .lines()
        .map(|l| l.with_context(|| format!("I/O error reading list file: {source}")))
        .filter_map(|l| match l {
            Err(e) => Some(Err(e)),
            Ok(line) => {
                let trimmed = line.trim().to_string();
                if trimmed.is_empty() || trimmed.starts_with('#') {
                    None
                } else {
                    Some(Ok(trimmed))
                }
            }
        })
        .map(|l| l.and_then(|l| parse_region_str(&l, flank)))
        .collect()
}

/// Parse a single region string: `chr:start-end` or `chr:pos+flank`.
pub fn parse_region_str(s: &str, flank: Option<u64>) -> Result<Region> {
    let (chr, rest) = s
        .split_once(':')
        .ok_or_else(|| anyhow!("Bad region (missing ':'): {}", s))?;

    // chr:start-end
    if let Some((start_s, end_s)) = rest.split_once('-') {
        let start: u64 = start_s
            .replace(',', "")
            .parse()
            .with_context(|| format!("Bad start in region: {s}"))?;
        let end: u64 = end_s
            .replace(',', "")
            .parse()
            .with_context(|| format!("Bad end in region: {s}"))?;
        // Coordinates are 1-based; reject zero values.
        if start == 0 {
            return Err(anyhow!(
                "start must be >= 1 (1-based coordinates) in region: {s}"
            ));
        }
        if end == 0 {
            return Err(anyhow!(
                "end must be >= 1 (1-based coordinates) in region: {s}"
            ));
        }
        return Ok(Region {
            name: None,
            chr: chr.to_string(),
            start: min(start, end),
            end: max(start, end),
        });
    }

    // chr:pos+flank
    let (pos_s, flank_s) = rest
        .split_once('+')
        .ok_or_else(|| anyhow!("Bad region format (expected start-end or pos+flank): {s}"))?;

    let pos: u64 = pos_s
        .replace(',', "")
        .parse()
        .with_context(|| format!("Bad position in region: {s}"))?;
    let inline_flank: u64 = flank_s
        .replace(',', "")
        .parse()
        .with_context(|| format!("Bad flank in region: {s}"))?;
    let effective_flank = inline_flank.max(flank.unwrap_or(0));

    Ok(Region {
        name: None,
        chr: chr.to_string(),
        start: pos.saturating_sub(effective_flank).max(1),
        end: pos.saturating_add(effective_flank),
    })
}

/// Parse a BED file (tab-separated: chr, start, end, optional name).
///
/// BED uses 0-based half-open coordinates. We convert to 1-based inclusive
/// by adding 1 to start (end stays the same, since half-open end equals
/// inclusive end in 1-based).
/// Skips comment lines (starting with #) and blank lines.
pub fn parse_regions_bed(path: &Path) -> Result<Vec<Region>> {
    let f =
        File::open(path).with_context(|| format!("Cannot open BED file: {}", path.display()))?;
    let reader = BufReader::new(f);
    let mut regions = Vec::new();

    for (line_num, line_result) in reader.lines().enumerate() {
        let line = line_result
            .with_context(|| format!("I/O error reading BED file: {}", path.display()))?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = trimmed.split('\t').collect();
        if fields.len() < 3 {
            return Err(anyhow!(
                "BED line {} has fewer than 3 fields: {}",
                line_num + 1,
                trimmed
            ));
        }
        let chr = fields[0].to_string();
        let start_0: u64 = fields[1]
            .parse()
            .with_context(|| format!("Bad start at BED line {}: {}", line_num + 1, fields[1]))?;
        let end: u64 = fields[2]
            .parse()
            .with_context(|| format!("Bad end at BED line {}: {}", line_num + 1, fields[2]))?;
        if end == 0 {
            return Err(anyhow!(
                "BED line {} has end=0, which produces an empty region",
                line_num + 1
            ));
        }
        let start_1based = start_0 + 1;
        if start_0 == end {
            return Err(anyhow!(
                "BED line {} has start == end ({}), which is an empty interval in \
                 0-based half-open coordinates",
                line_num + 1,
                start_0
            ));
        }
        if start_1based > end {
            return Err(anyhow!(
                "BED line {} is malformed: after converting to 1-based coordinates, \
                 start ({}) > end ({})",
                line_num + 1,
                start_1based,
                end
            ));
        }
        let name = if fields.len() >= 4 && !fields[3].is_empty() {
            Some(fields[3].to_string())
        } else {
            None
        };
        regions.push(Region {
            name,
            chr,
            start: start_1based,
            end,
        });
    }
    Ok(regions)
}

/// Remove duplicate regions (same chr, start, end). Returns the number of duplicates removed.
/// Preserves the first occurrence of each unique region.
#[allow(dead_code)]
pub fn deduplicate_regions(regions: &mut Vec<Region>) -> usize {
    use std::collections::HashSet;
    let original_len = regions.len();
    let mut seen = HashSet::new();
    regions.retain(|r| seen.insert((r.chr.clone(), r.start, r.end)));
    original_len - regions.len()
}

/// Sort regions by chromosome (natural order) then start position.
/// Natural order means chr1, chr2, ..., chr10 rather than chr1, chr10, chr2.
#[allow(dead_code)]
pub fn sort_regions(regions: &mut [Region]) {
    regions.sort_by(|a, b| natural_chr_cmp(&a.chr, &b.chr).then(a.start.cmp(&b.start)));
}

/// Compare two chromosome names using natural sort order.
/// Splits each name into text and numeric segments and compares them pairwise.
fn natural_chr_cmp(a: &str, b: &str) -> std::cmp::Ordering {
    let mut ai = a.chars().peekable();
    let mut bi = b.chars().peekable();

    loop {
        match (ai.peek(), bi.peek()) {
            (None, None) => return std::cmp::Ordering::Equal,
            (None, Some(_)) => return std::cmp::Ordering::Less,
            (Some(_), None) => return std::cmp::Ordering::Greater,
            (Some(&ac), Some(&bc)) => {
                if ac.is_ascii_digit() && bc.is_ascii_digit() {
                    // Compare numeric segments.
                    let an = consume_number(&mut ai);
                    let bn = consume_number(&mut bi);
                    match an.cmp(&bn) {
                        std::cmp::Ordering::Equal => continue,
                        other => return other,
                    }
                } else {
                    ai.next();
                    bi.next();
                    match ac.cmp(&bc) {
                        std::cmp::Ordering::Equal => continue,
                        other => return other,
                    }
                }
            }
        }
    }
}

/// Consume consecutive digits from a peekable char iterator and return as u64.
fn consume_number(it: &mut std::iter::Peekable<std::str::Chars>) -> u64 {
    let mut n: u64 = 0;
    while let Some(&c) = it.peek() {
        if c.is_ascii_digit() {
            n = n.saturating_mul(10).saturating_add(c as u64 - '0' as u64);
            it.next();
        } else {
            break;
        }
    }
    n
}

// ─── Unit tests ──────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── parse_region_str ─────────────────────────────────────────────────

    #[test]
    fn parse_region_str_range() {
        let r = parse_region_str("chr1:100-200", None).unwrap();
        assert_eq!(r.chr, "chr1");
        assert_eq!(r.start, 100);
        assert_eq!(r.end, 200);
        assert!(r.name.is_none());
    }

    #[test]
    fn parse_region_str_swapped() {
        let r = parse_region_str("chr1:200-100", None).unwrap();
        assert_eq!(r.start, 100);
        assert_eq!(r.end, 200);
    }

    #[test]
    fn parse_region_str_with_commas() {
        let r = parse_region_str("chr1:1,000-2,000", None).unwrap();
        assert_eq!(r.start, 1000);
        assert_eq!(r.end, 2000);
    }

    #[test]
    fn parse_region_str_position_flank_inline() {
        let r = parse_region_str("chr1:1000+500", None).unwrap();
        assert_eq!(r.chr, "chr1");
        assert_eq!(r.start, 500);
        assert_eq!(r.end, 1500);
    }

    #[test]
    fn parse_region_str_position_flank_cli_override() {
        // CLI flank is larger than inline flank, so it wins.
        let r = parse_region_str("chr1:1000+100", Some(500)).unwrap();
        assert_eq!(r.start, 500);
        assert_eq!(r.end, 1500);
    }

    #[test]
    fn parse_region_str_position_clamps_to_one() {
        let r = parse_region_str("chr1:3+10", None).unwrap();
        assert_eq!(r.start, 1);
        assert_eq!(r.end, 13);
    }

    #[test]
    fn parse_region_str_missing_colon() {
        assert!(parse_region_str("chr1_100_200", None).is_err());
    }

    #[test]
    fn parse_region_str_bad_start() {
        assert!(parse_region_str("chr1:abc-200", None).is_err());
    }

    #[test]
    fn parse_region_str_bad_format() {
        assert!(parse_region_str("chr1:100", None).is_err());
    }

    // ── parse_regions_inline ─────────────────────────────────────────────

    #[test]
    fn parse_regions_inline_multiple() {
        let regions = parse_regions_inline("chr1:1-10,chr2:20-30", None).unwrap();
        assert_eq!(regions.len(), 2);
        assert_eq!(regions[0].chr, "chr1");
        assert_eq!(regions[1].chr, "chr2");
    }

    #[test]
    fn parse_regions_inline_trailing_comma() {
        let regions = parse_regions_inline("chr1:1-10,", None).unwrap();
        assert_eq!(regions.len(), 1);
    }

    #[test]
    fn parse_regions_inline_empty() {
        let regions = parse_regions_inline("", None).unwrap();
        assert_eq!(regions.len(), 0);
    }

    // ── overflow / saturation ──────────────────────────────────────────

    #[test]
    fn parse_region_str_large_position_flank_saturates() {
        // A position near u64::MAX with a flank must saturate to u64::MAX
        // rather than wrapping around to a small value.
        let r = parse_region_str(&format!("chr1:{}+100", u64::MAX - 10), None).unwrap();
        assert_eq!(r.end, u64::MAX, "end should saturate to u64::MAX, not wrap");
        // start should also saturate sensibly (stay above 1)
        assert!(r.start >= 1);
    }

    // ── zero-coordinate rejection ─────────────────────────────────────────

    #[test]
    fn parse_region_str_zero_start_errors() {
        let err = parse_region_str("chr1:0-10", None).unwrap_err();
        let msg = err.to_string();
        assert!(
            msg.contains("1-based coordinates"),
            "expected 1-based coordinates message, got: {msg}"
        );
    }

    #[test]
    fn parse_region_str_zero_end_errors() {
        let err = parse_region_str("chr1:10-0", None).unwrap_err();
        let msg = err.to_string();
        assert!(
            msg.contains("1-based coordinates"),
            "expected 1-based coordinates message, got: {msg}"
        );
    }

    // ── deduplicate_regions ─────────────────────────────────────────────

    #[test]
    fn dedup_removes_exact_duplicates() {
        let mut regions = vec![
            Region {
                name: None,
                chr: "chr1".into(),
                start: 1,
                end: 10,
            },
            Region {
                name: Some("foo".into()),
                chr: "chr1".into(),
                start: 1,
                end: 10,
            },
            Region {
                name: None,
                chr: "chr2".into(),
                start: 1,
                end: 10,
            },
        ];
        let removed = deduplicate_regions(&mut regions);
        assert_eq!(removed, 1);
        assert_eq!(regions.len(), 2);
        assert_eq!(regions[0].chr, "chr1");
        assert_eq!(regions[1].chr, "chr2");
    }

    #[test]
    fn dedup_no_duplicates() {
        let mut regions = vec![
            Region {
                name: None,
                chr: "chr1".into(),
                start: 1,
                end: 10,
            },
            Region {
                name: None,
                chr: "chr2".into(),
                start: 1,
                end: 10,
            },
        ];
        let removed = deduplicate_regions(&mut regions);
        assert_eq!(removed, 0);
        assert_eq!(regions.len(), 2);
    }

    // ── sort_regions ────────────────────────────────────────────────────

    // ── parse_regions_bed ──────────────────────────────────────────────

    #[test]
    fn bed_valid_region() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        std::io::Write::write_all(&mut tmp.as_file(), b"chr1\t100\t200\n").unwrap();
        let regions = parse_regions_bed(tmp.path()).unwrap();
        assert_eq!(regions.len(), 1);
        // BED 0-based half-open [100, 200) -> 1-based inclusive [101, 200]
        assert_eq!(regions[0].start, 101);
        assert_eq!(regions[0].end, 200);
    }

    #[test]
    fn bed_start_equals_end_errors() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        // start == end == 5 in BED means an empty interval.
        std::io::Write::write_all(&mut tmp.as_file(), b"chr1\t5\t5\n").unwrap();
        let err = parse_regions_bed(tmp.path()).unwrap_err();
        let msg = err.to_string();
        assert!(
            msg.contains("empty interval"),
            "expected empty interval message, got: {msg}"
        );
    }

    #[test]
    fn bed_start_greater_than_end_errors() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        // start=10, end=5 in BED is malformed.
        std::io::Write::write_all(&mut tmp.as_file(), b"chr1\t10\t5\n").unwrap();
        let err = parse_regions_bed(tmp.path()).unwrap_err();
        let msg = err.to_string();
        assert!(
            msg.contains("malformed"),
            "expected malformed message, got: {msg}"
        );
    }

    // ── sort_regions ────────────────────────────────────────────────────

    #[test]
    fn sort_natural_chromosome_order() {
        let mut regions = vec![
            Region {
                name: None,
                chr: "chr10".into(),
                start: 1,
                end: 10,
            },
            Region {
                name: None,
                chr: "chr2".into(),
                start: 1,
                end: 10,
            },
            Region {
                name: None,
                chr: "chr1".into(),
                start: 20,
                end: 30,
            },
            Region {
                name: None,
                chr: "chr1".into(),
                start: 1,
                end: 10,
            },
        ];
        sort_regions(&mut regions);
        assert_eq!(regions[0].chr, "chr1");
        assert_eq!(regions[0].start, 1);
        assert_eq!(regions[1].chr, "chr1");
        assert_eq!(regions[1].start, 20);
        assert_eq!(regions[2].chr, "chr2");
        assert_eq!(regions[3].chr, "chr10");
    }
}
