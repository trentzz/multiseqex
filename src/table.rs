use anyhow::{Context, Result, anyhow};
use std::cmp::{max, min};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use crate::region::{Region, resolve_flanks};

/// Build a case-insensitive header-name to column-index map.
pub fn build_header_map(headers: &csv::StringRecord) -> HashMap<String, usize> {
    headers
        .iter()
        .enumerate()
        .map(|(i, h)| (h.trim().to_uppercase(), i))
        .collect()
}

/// Parse a user-supplied delimiter string into a byte.
///
/// Accepts "tab", "comma", or a single ASCII character.
pub fn parse_delimiter_flag(s: &str) -> Result<u8> {
    match s.to_ascii_lowercase().as_str() {
        "tab" => Ok(b'\t'),
        "comma" => Ok(b','),
        _ => {
            let bytes = s.as_bytes();
            if bytes.len() == 1 {
                Ok(bytes[0])
            } else {
                Err(anyhow!(
                    "Invalid --delimiter value \"{s}\": expected \"tab\", \"comma\", or a single character"
                ))
            }
        }
    }
}

/// Detect the delimiter for a table file.
///
/// Priority:
/// 1. Explicit override from `--delimiter`.
/// 2. File extension: `.tsv` leads to tab.
/// 3. Content sniffing: if the first line contains a tab, use tab.
/// 4. Default to comma.
pub fn detect_delimiter(path: &Path, cli_delimiter: Option<&str>) -> Result<u8> {
    // 1. Explicit override.
    if let Some(s) = cli_delimiter {
        return parse_delimiter_flag(s);
    }

    // 2. Extension check.
    let is_tsv = path
        .extension()
        .and_then(|e| e.to_str())
        .is_some_and(|e| e.eq_ignore_ascii_case("tsv"));
    if is_tsv {
        return Ok(b'\t');
    }

    // 3. Sniff the first line for tabs.
    if let Ok(file) = File::open(path) {
        let reader = BufReader::new(file);
        if let Some(Ok(first_line)) = reader.lines().next() {
            if first_line.contains('\t') {
                return Ok(b'\t');
            }
        }
    }

    // 4. Default to comma.
    Ok(b',')
}

/// Parse a numeric field from a CSV record, with a contextual error message.
fn parse_u64_field(
    rec: &csv::StringRecord,
    col_idx: usize,
    field_name: &str,
    row: usize,
) -> Result<u64> {
    rec.get(col_idx)
        .ok_or_else(|| anyhow!("Missing {field_name} at row {row}"))?
        .trim()
        .replace(',', "")
        .parse()
        .with_context(|| format!("Bad {field_name} at row {row}"))
}

/// Look up a required column by name, returning a helpful error if absent.
fn require_column(
    hmap: &HashMap<String, usize>,
    name: &str,
    headers: &csv::StringRecord,
) -> Result<usize> {
    hmap.get(name)
        .copied()
        .ok_or_else(|| anyhow!("Table missing required {name} column (found: {headers:?})"))
}

/// Read a string field from a CSV record.
fn read_string_field(rec: &csv::StringRecord, col_idx: usize) -> Result<String> {
    Ok(rec.get(col_idx).unwrap_or("").trim().to_string())
}

/// Read an optional NAME field (returns `None` if the column is absent or empty).
fn read_optional_name(rec: &csv::StringRecord, name_idx: Option<usize>) -> Option<String> {
    name_idx.and_then(|idx| {
        rec.get(idx)
            .map(|v| v.trim().to_string())
            .filter(|v| !v.is_empty())
    })
}

enum TableMode {
    Range { start_idx: usize, end_idx: usize },
    Position { pos_idx: usize },
}

/// Read an optional STRAND field (returns `None` if the column is absent or not +/-/.).
fn read_optional_strand(rec: &csv::StringRecord, strand_idx: Option<usize>) -> Option<char> {
    strand_idx.and_then(|idx| {
        rec.get(idx).and_then(|v| match v.trim() {
            "+" => Some('+'),
            "-" => Some('-'),
            "." => Some('.'),
            _ => None,
        })
    })
}

/// Parse a CSV/TSV table with named columns: CHROM, START, END, POS, NAME, STRAND.
pub fn parse_regions_table(
    path: &Path,
    flank: Option<u64>,
    flank_left: Option<u64>,
    flank_right: Option<u64>,
    cli_delimiter: Option<&str>,
) -> Result<Vec<Region>> {
    let delim = detect_delimiter(path, cli_delimiter)?;
    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(true)
        .delimiter(delim)
        .from_path(path)
        .with_context(|| format!("Cannot open table: {}", path.display()))?;

    let headers = rdr.headers()?.clone();
    let hmap = build_header_map(&headers);

    let chrom_idx = require_column(&hmap, "CHROM", &headers)?;
    let name_idx = hmap.get("NAME").copied();
    let strand_idx = hmap.get("STRAND").copied();

    let has_start = hmap.get("START").copied();
    let has_end = hmap.get("END").copied();
    let has_pos = hmap.get("POS").copied();

    let mode = match (has_start, has_end, has_pos) {
        (Some(s), Some(e), None) => TableMode::Range {
            start_idx: s,
            end_idx: e,
        },
        (None, None, Some(p)) => TableMode::Position { pos_idx: p },
        (Some(_), Some(_), Some(_)) => {
            return Err(anyhow!(
                "Table has both START/END and POS columns — ambiguous. \
                 Use either START+END (range) or POS (position)."
            ));
        }
        (Some(_), None, _) | (None, Some(_), _) => {
            return Err(anyhow!(
                "Table has only one of START/END — both are required for range mode (found: {headers:?})"
            ));
        }
        _ => {
            return Err(anyhow!(
                "Table must have START+END or POS columns (found: {headers:?})"
            ));
        }
    };

    if matches!(mode, TableMode::Position { .. }) && flank.is_none() && flank_left.is_none() {
        return Err(anyhow!(
            "--flank is required when table uses POS column (position mode)"
        ));
    }
    let (fl, fr) = resolve_flanks(flank, flank_left, flank_right);

    let mut out = Vec::new();
    for (i, rec) in rdr.records().enumerate() {
        let rec = rec?;
        let row = i + 2;
        let chr = read_string_field(&rec, chrom_idx)?;
        let name = read_optional_name(&rec, name_idx);
        let strand = read_optional_strand(&rec, strand_idx);

        let region = match &mode {
            TableMode::Range { start_idx, end_idx } => {
                let s = parse_u64_field(&rec, *start_idx, "START", row)?;
                let e = parse_u64_field(&rec, *end_idx, "END", row)?;
                // Coordinates are 1-based; reject zero values.
                if s == 0 {
                    return Err(anyhow!(
                        "START must be >= 1 (1-based coordinates) at row {row}"
                    ));
                }
                if e == 0 {
                    return Err(anyhow!(
                        "END must be >= 1 (1-based coordinates) at row {row}"
                    ));
                }
                Region {
                    name,
                    chr,
                    start: min(s, e),
                    end: max(s, e),
                    strand,
                }
            }
            TableMode::Position { pos_idx } => {
                let p = parse_u64_field(&rec, *pos_idx, "POS", row)?;
                // Coordinates are 1-based; reject zero values.
                if p == 0 {
                    return Err(anyhow!(
                        "POS must be >= 1 (1-based coordinates) at row {row}"
                    ));
                }
                Region {
                    name,
                    chr,
                    start: p.saturating_sub(fl).max(1),
                    end: p.saturating_add(fr),
                    strand,
                }
            }
        };
        out.push(region);
    }
    Ok(out)
}

enum SvMode {
    Range {
        start_left_idx: usize,
        end_left_idx: usize,
        start_right_idx: usize,
        end_right_idx: usize,
    },
    Position {
        pos_left_idx: usize,
        pos_right_idx: usize,
    },
}

/// Parse a CSV/TSV SV table with named columns.
pub fn parse_regions_sv_table(
    path: &Path,
    flank: Option<u64>,
    flank_left: Option<u64>,
    flank_right: Option<u64>,
    cli_delimiter: Option<&str>,
) -> Result<Vec<Region>> {
    let delim = detect_delimiter(path, cli_delimiter)?;
    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(true)
        .delimiter(delim)
        .from_path(path)
        .with_context(|| format!("Cannot open SV table: {}", path.display()))?;

    let headers = rdr.headers()?.clone();
    let hmap = build_header_map(&headers);

    let chrom_left_idx = require_column(&hmap, "CHROM_LEFT", &headers)?;
    let chrom_right_idx = require_column(&hmap, "CHROM_RIGHT", &headers)?;
    let name_idx = hmap.get("NAME").copied();
    let strand_idx = hmap.get("STRAND").copied();

    let has_sl = hmap.get("START_LEFT").copied();
    let has_el = hmap.get("END_LEFT").copied();
    let has_sr = hmap.get("START_RIGHT").copied();
    let has_er = hmap.get("END_RIGHT").copied();
    let has_pl = hmap.get("POS_LEFT").copied();
    let has_pr = hmap.get("POS_RIGHT").copied();

    let mode = match (has_sl, has_el, has_sr, has_er, has_pl, has_pr) {
        (Some(sl), Some(el), Some(sr), Some(er), _, _) => SvMode::Range {
            start_left_idx: sl,
            end_left_idx: el,
            start_right_idx: sr,
            end_right_idx: er,
        },
        (None, None, None, None, Some(pl), Some(pr)) => SvMode::Position {
            pos_left_idx: pl,
            pos_right_idx: pr,
        },
        _ => {
            return Err(anyhow!(
                "SV table must have START_LEFT+END_LEFT+START_RIGHT+END_RIGHT (range) \
                 or POS_LEFT+POS_RIGHT (position). Found: {headers:?}"
            ));
        }
    };

    if matches!(mode, SvMode::Position { .. }) && flank.is_none() && flank_left.is_none() {
        return Err(anyhow!(
            "--flank is required when SV table uses POS_LEFT/POS_RIGHT (position mode)"
        ));
    }
    let (fl, fr) = resolve_flanks(flank, flank_left, flank_right);

    let mut out = Vec::new();
    for (i, rec) in rdr.records().enumerate() {
        let rec = rec?;
        let row = i + 2;
        let chr_left = read_string_field(&rec, chrom_left_idx)?;
        let chr_right = read_string_field(&rec, chrom_right_idx)?;
        let name = read_optional_name(&rec, name_idx);
        let strand = read_optional_strand(&rec, strand_idx);

        match &mode {
            SvMode::Range {
                start_left_idx,
                end_left_idx,
                start_right_idx,
                end_right_idx,
            } => {
                let sl = parse_u64_field(&rec, *start_left_idx, "START_LEFT", row)?;
                let el = parse_u64_field(&rec, *end_left_idx, "END_LEFT", row)?;
                let sr = parse_u64_field(&rec, *start_right_idx, "START_RIGHT", row)?;
                let er = parse_u64_field(&rec, *end_right_idx, "END_RIGHT", row)?;
                // Coordinates are 1-based; reject zero values.
                for (val, label) in [
                    (sl, "START_LEFT"),
                    (el, "END_LEFT"),
                    (sr, "START_RIGHT"),
                    (er, "END_RIGHT"),
                ] {
                    if val == 0 {
                        return Err(anyhow!(
                            "{label} must be >= 1 (1-based coordinates) at row {row}"
                        ));
                    }
                }
                out.push(Region {
                    name: name.clone(),
                    chr: chr_left,
                    start: min(sl, el),
                    end: max(sl, el),
                    strand,
                });
                out.push(Region {
                    name,
                    chr: chr_right,
                    start: min(sr, er),
                    end: max(sr, er),
                    strand,
                });
            }
            SvMode::Position {
                pos_left_idx,
                pos_right_idx,
            } => {
                let pl = parse_u64_field(&rec, *pos_left_idx, "POS_LEFT", row)?;
                let pr = parse_u64_field(&rec, *pos_right_idx, "POS_RIGHT", row)?;
                // Coordinates are 1-based; reject zero values.
                for (val, label) in [(pl, "POS_LEFT"), (pr, "POS_RIGHT")] {
                    if val == 0 {
                        return Err(anyhow!(
                            "{label} must be >= 1 (1-based coordinates) at row {row}"
                        ));
                    }
                }
                out.push(Region {
                    name: name.clone(),
                    chr: chr_left,
                    start: pl.saturating_sub(fl).max(1),
                    end: pl.saturating_add(fr),
                    strand,
                });
                out.push(Region {
                    name,
                    chr: chr_right,
                    start: pr.saturating_sub(fl).max(1),
                    end: pr.saturating_add(fr),
                    strand,
                });
            }
        }
    }
    Ok(out)
}

// ─── Unit tests ──────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── detect_delimiter ─────────────────────────────────────────────────

    #[test]
    fn detect_delimiter_csv() {
        assert_eq!(detect_delimiter(Path::new("file.csv"), None).unwrap(), b',');
    }

    #[test]
    fn detect_delimiter_tsv() {
        assert_eq!(
            detect_delimiter(Path::new("file.tsv"), None).unwrap(),
            b'\t'
        );
    }

    #[test]
    fn detect_delimiter_tsv_uppercase() {
        assert_eq!(
            detect_delimiter(Path::new("file.TSV"), None).unwrap(),
            b'\t'
        );
    }

    #[test]
    fn detect_delimiter_no_extension() {
        assert_eq!(detect_delimiter(Path::new("file"), None).unwrap(), b',');
    }

    #[test]
    fn detect_delimiter_sniff_tabs_in_txt() {
        // A .txt file whose first line contains tabs should be detected as tab-delimited.
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("data.txt");
        std::fs::write(&path, "CHROM\tSTART\tEND\nchr1\t1\t10\n").unwrap();
        assert_eq!(detect_delimiter(&path, None).unwrap(), b'\t');
    }

    #[test]
    fn detect_delimiter_sniff_comma_in_txt() {
        // A .txt file with no tabs defaults to comma.
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("data.txt");
        std::fs::write(&path, "CHROM,START,END\nchr1,1,10\n").unwrap();
        assert_eq!(detect_delimiter(&path, None).unwrap(), b',');
    }

    #[test]
    fn detect_delimiter_cli_override_tab() {
        assert_eq!(
            detect_delimiter(Path::new("file.csv"), Some("tab")).unwrap(),
            b'\t'
        );
    }

    #[test]
    fn detect_delimiter_cli_override_comma() {
        assert_eq!(
            detect_delimiter(Path::new("file.tsv"), Some("comma")).unwrap(),
            b','
        );
    }

    #[test]
    fn detect_delimiter_cli_override_single_char() {
        assert_eq!(
            detect_delimiter(Path::new("file.csv"), Some(";")).unwrap(),
            b';'
        );
    }

    #[test]
    fn detect_delimiter_cli_override_invalid() {
        assert!(detect_delimiter(Path::new("file.csv"), Some("abc")).is_err());
    }

    // ── build_header_map ─────────────────────────────────────────────────

    #[test]
    fn build_header_map_case_insensitive() {
        let rec = csv::StringRecord::from(vec!["chrom", "Start", "END"]);
        let map = build_header_map(&rec);
        assert!(map.contains_key("CHROM"));
        assert!(map.contains_key("START"));
        assert!(map.contains_key("END"));
    }

    #[test]
    fn build_header_map_trims_whitespace() {
        let rec = csv::StringRecord::from(vec![" CHROM ", " START"]);
        let map = build_header_map(&rec);
        assert!(map.contains_key("CHROM"));
        assert!(map.contains_key("START"));
    }
}
