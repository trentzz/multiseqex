use anyhow::{Context, Result, anyhow};
use std::cmp::{max, min};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// A genomic interval (1-based, inclusive on both ends).
#[derive(Debug, Clone)]
pub(crate) struct Region {
    pub(crate) name: Option<String>,
    pub(crate) chr: String,
    pub(crate) start: u64,
    pub(crate) end: u64,
}

/// Parse comma-separated inline region strings.
pub(crate) fn parse_regions_inline(s: &str, flank: Option<u64>) -> Result<Vec<Region>> {
    s.split(',')
        .filter(|t| !t.trim().is_empty())
        .map(|t| parse_region_str(t.trim(), flank))
        .collect()
}

/// Parse a file with one region string per line.
pub(crate) fn parse_regions_list(path: &Path, flank: Option<u64>) -> Result<Vec<Region>> {
    let f =
        File::open(path).with_context(|| format!("Cannot open list file: {}", path.display()))?;
    BufReader::new(f)
        .lines()
        .map(|l| l.with_context(|| format!("I/O error reading list file: {}", path.display())))
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
pub(crate) fn parse_region_str(s: &str, flank: Option<u64>) -> Result<Region> {
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
}
