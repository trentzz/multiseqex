//! Interval operations on genomic regions.
//!
//! Provides `subtract` and `intersect` operations that modify a list of
//! regions based on intervals from a BED file. Both operate on regions
//! BEFORE extraction: they modify the region list, potentially splitting
//! or trimming regions.

use anyhow::{Context, Result};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use crate::region::Region;

/// A simple half-open interval [start, end) used internally for arithmetic.
/// Stored as 1-based inclusive to match `Region` coordinates.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Interval {
    start: u64,
    end: u64,
}

/// Load a BED file into per-chromosome sorted interval lists.
/// BED coordinates are 0-based half-open; we convert to 1-based inclusive
/// to match `Region`.
fn load_bed_intervals(path: &Path) -> Result<HashMap<String, Vec<Interval>>> {
    let f =
        File::open(path).with_context(|| format!("Cannot open BED file: {}", path.display()))?;
    let reader = BufReader::new(f);
    let mut map: HashMap<String, Vec<Interval>> = HashMap::new();

    for (line_num, line_result) in reader.lines().enumerate() {
        let line =
            line_result.with_context(|| format!("I/O error reading BED: {}", path.display()))?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = trimmed.split('\t').collect();
        if fields.len() < 3 {
            return Err(anyhow::anyhow!(
                "BED line {} has fewer than 3 fields: {}",
                line_num + 1,
                trimmed
            ));
        }
        let chr = fields[0].to_string();
        let start_0: u64 = fields[1]
            .parse()
            .with_context(|| format!("Bad start at BED line {}: {}", line_num + 1, fields[1]))?;
        let end_0: u64 = fields[2]
            .parse()
            .with_context(|| format!("Bad end at BED line {}: {}", line_num + 1, fields[2]))?;
        // Convert 0-based half-open to 1-based inclusive.
        let start = start_0 + 1;
        let end = end_0;
        if start > end {
            continue; // Skip empty intervals.
        }
        map.entry(chr).or_default().push(Interval { start, end });
    }

    // Sort intervals per chromosome by start position.
    for intervals in map.values_mut() {
        intervals.sort_by_key(|iv| (iv.start, iv.end));
    }

    Ok(map)
}

/// Subtract intervals from a set of regions.
///
/// For each input region, removes any portions that overlap with intervals
/// from the BED file. A single region may be split into multiple pieces.
pub fn subtract_regions(regions: &[Region], bed_path: &Path) -> Result<Vec<Region>> {
    let bed_intervals = load_bed_intervals(bed_path)?;
    let mut result = Vec::new();

    for r in regions {
        let Some(intervals) = bed_intervals.get(&r.chr) else {
            // No intervals for this chromosome: keep region unchanged.
            result.push(r.clone());
            continue;
        };

        // Find overlapping intervals and subtract them.
        let pieces = subtract_from_interval(r.start, r.end, intervals);
        for (i, (s, e)) in pieces.iter().enumerate() {
            let mut region = r.clone();
            region.start = *s;
            region.end = *e;
            if pieces.len() > 1 {
                // Append piece index to name when a region is split.
                let base = r.name.as_deref().unwrap_or("");
                let suffix = format!(
                    "{}_{}_{}",
                    if base.is_empty() { "" } else { "_" },
                    "part",
                    i + 1
                );
                region.name = Some(format!(
                    "{}{}",
                    r.name
                        .as_deref()
                        .unwrap_or(&format!("{}:{}-{}", r.chr, r.start, r.end)),
                    suffix
                ));
            }
            result.push(region);
        }
    }

    Ok(result)
}

/// Subtract sorted intervals from a single [start, end] range.
/// Returns the remaining pieces.
fn subtract_from_interval(start: u64, end: u64, intervals: &[Interval]) -> Vec<(u64, u64)> {
    let mut pieces = Vec::new();
    let mut cursor = start;

    for iv in intervals {
        // Skip intervals entirely before our range.
        if iv.end < cursor {
            continue;
        }
        // Stop if the interval starts after our range.
        if iv.start > end {
            break;
        }
        // There is overlap. Keep the portion before the interval.
        if cursor < iv.start {
            pieces.push((cursor, iv.start - 1));
        }
        // Move cursor past the interval.
        cursor = iv.end + 1;
    }

    // Keep the remaining portion after the last interval.
    if cursor <= end {
        pieces.push((cursor, end));
    }

    pieces
}

/// Intersect regions with intervals from a BED file.
///
/// For each input region, keeps only the portions that overlap with
/// intervals in the BED file. A single region may produce multiple
/// output regions.
pub fn intersect_regions(regions: &[Region], bed_path: &Path) -> Result<Vec<Region>> {
    let bed_intervals = load_bed_intervals(bed_path)?;
    let mut result = Vec::new();

    for r in regions {
        let Some(intervals) = bed_intervals.get(&r.chr) else {
            // No intervals for this chromosome: region is removed entirely.
            continue;
        };

        let pieces = intersect_with_intervals(r.start, r.end, intervals);
        for (i, (s, e)) in pieces.iter().enumerate() {
            let mut region = r.clone();
            region.start = *s;
            region.end = *e;
            if pieces.len() > 1 {
                let base = r.name.as_deref().unwrap_or("");
                let suffix = format!(
                    "{}_{}_{}",
                    if base.is_empty() { "" } else { "_" },
                    "part",
                    i + 1
                );
                region.name = Some(format!(
                    "{}{}",
                    r.name
                        .as_deref()
                        .unwrap_or(&format!("{}:{}-{}", r.chr, r.start, r.end)),
                    suffix
                ));
            }
            result.push(region);
        }
    }

    Ok(result)
}

/// Intersect a single [start, end] range with sorted intervals.
/// Returns the overlapping portions.
fn intersect_with_intervals(start: u64, end: u64, intervals: &[Interval]) -> Vec<(u64, u64)> {
    let mut pieces = Vec::new();

    for iv in intervals {
        // Skip intervals entirely before our range.
        if iv.end < start {
            continue;
        }
        // Stop if the interval starts after our range.
        if iv.start > end {
            break;
        }
        // Compute the overlap.
        let overlap_start = start.max(iv.start);
        let overlap_end = end.min(iv.end);
        if overlap_start <= overlap_end {
            pieces.push((overlap_start, overlap_end));
        }
    }

    pieces
}

// ─── Unit tests ──────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn make_region(chr: &str, start: u64, end: u64, name: Option<&str>) -> Region {
        Region {
            name: name.map(|s| s.to_string()),
            chr: chr.to_string(),
            start,
            end,
            strand: None,
            alt_info: None,
        }
    }

    // ── subtract_from_interval ───────────────────────────────────────────

    #[test]
    fn subtract_no_overlap() {
        let intervals = vec![Interval { start: 50, end: 60 }];
        let pieces = subtract_from_interval(1, 10, &intervals);
        assert_eq!(pieces, vec![(1, 10)]);
    }

    #[test]
    fn subtract_complete_overlap() {
        let intervals = vec![Interval { start: 1, end: 100 }];
        let pieces = subtract_from_interval(10, 50, &intervals);
        assert!(pieces.is_empty());
    }

    #[test]
    fn subtract_splits_into_two() {
        let intervals = vec![Interval { start: 20, end: 30 }];
        let pieces = subtract_from_interval(10, 50, &intervals);
        assert_eq!(pieces, vec![(10, 19), (31, 50)]);
    }

    #[test]
    fn subtract_trims_start() {
        let intervals = vec![Interval { start: 1, end: 20 }];
        let pieces = subtract_from_interval(10, 50, &intervals);
        assert_eq!(pieces, vec![(21, 50)]);
    }

    #[test]
    fn subtract_trims_end() {
        let intervals = vec![Interval { start: 40, end: 60 }];
        let pieces = subtract_from_interval(10, 50, &intervals);
        assert_eq!(pieces, vec![(10, 39)]);
    }

    #[test]
    fn subtract_multiple_intervals() {
        let intervals = vec![
            Interval { start: 15, end: 20 },
            Interval { start: 30, end: 35 },
        ];
        let pieces = subtract_from_interval(10, 50, &intervals);
        assert_eq!(pieces, vec![(10, 14), (21, 29), (36, 50)]);
    }

    // ── intersect_with_intervals ─────────────────────────────────────────

    #[test]
    fn intersect_no_overlap() {
        let intervals = vec![Interval { start: 50, end: 60 }];
        let pieces = intersect_with_intervals(1, 10, &intervals);
        assert!(pieces.is_empty());
    }

    #[test]
    fn intersect_complete_overlap() {
        let intervals = vec![Interval { start: 1, end: 100 }];
        let pieces = intersect_with_intervals(10, 50, &intervals);
        assert_eq!(pieces, vec![(10, 50)]);
    }

    #[test]
    fn intersect_partial_overlap() {
        let intervals = vec![Interval { start: 5, end: 15 }];
        let pieces = intersect_with_intervals(10, 50, &intervals);
        assert_eq!(pieces, vec![(10, 15)]);
    }

    #[test]
    fn intersect_multiple_intervals() {
        let intervals = vec![
            Interval { start: 5, end: 15 },
            Interval { start: 30, end: 40 },
            Interval { start: 60, end: 70 },
        ];
        let pieces = intersect_with_intervals(10, 50, &intervals);
        assert_eq!(pieces, vec![(10, 15), (30, 40)]);
    }

    // ── subtract_regions (file-based) ────────────────────────────────────

    #[test]
    fn subtract_regions_from_bed() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        std::io::Write::write_all(
            &mut tmp.as_file(),
            b"chr1\t19\t30\n", // 0-based [19,30) -> 1-based [20,30]
        )
        .unwrap();

        let regions = vec![make_region("chr1", 10, 50, None)];
        let result = subtract_regions(&regions, tmp.path()).unwrap();
        assert_eq!(result.len(), 2);
        assert_eq!(result[0].start, 10);
        assert_eq!(result[0].end, 19);
        assert_eq!(result[1].start, 31);
        assert_eq!(result[1].end, 50);
    }

    #[test]
    fn subtract_regions_no_overlap() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        std::io::Write::write_all(&mut tmp.as_file(), b"chr2\t0\t100\n").unwrap();

        let regions = vec![make_region("chr1", 10, 50, None)];
        let result = subtract_regions(&regions, tmp.path()).unwrap();
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].start, 10);
        assert_eq!(result[0].end, 50);
    }

    // ── intersect_regions (file-based) ───────────────────────────────────

    #[test]
    fn intersect_regions_from_bed() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        std::io::Write::write_all(
            &mut tmp.as_file(),
            b"chr1\t14\t25\n", // 0-based [14,25) -> 1-based [15,25]
        )
        .unwrap();

        let regions = vec![make_region("chr1", 10, 50, None)];
        let result = intersect_regions(&regions, tmp.path()).unwrap();
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].start, 15);
        assert_eq!(result[0].end, 25);
    }

    #[test]
    fn intersect_regions_no_overlap() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        std::io::Write::write_all(&mut tmp.as_file(), b"chr2\t0\t100\n").unwrap();

        let regions = vec![make_region("chr1", 10, 50, None)];
        let result = intersect_regions(&regions, tmp.path()).unwrap();
        assert!(result.is_empty());
    }
}
