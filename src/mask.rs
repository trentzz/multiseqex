//! Sequence masking using BED intervals.
//!
//! After extraction, overlapping intervals from a mask BED file are applied
//! to replace or lowercase the affected bases.

use anyhow::{Context, Result, anyhow};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// Masking mode.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum MaskMode {
    /// Replace masked bases with 'N'.
    #[default]
    Hard,
    /// Lowercase masked bases.
    Soft,
}

/// A collection of mask intervals keyed by chromosome.
///
/// Intervals are stored as 1-based inclusive `(start, end)` pairs, sorted
/// by start position within each chromosome for efficient overlap scanning.
#[derive(Debug, Clone)]
pub struct MaskIndex {
    intervals: HashMap<String, Vec<(u64, u64)>>,
}

impl MaskIndex {
    /// Load a mask BED file. BED coordinates are 0-based half-open and are
    /// converted to 1-based inclusive internally (same as region coordinates).
    pub fn from_bed(path: &Path) -> Result<Self> {
        let f = File::open(path)
            .with_context(|| format!("Cannot open mask BED file: {}", path.display()))?;
        let reader = BufReader::new(f);
        let mut intervals: HashMap<String, Vec<(u64, u64)>> = HashMap::new();

        for (line_num, line_result) in reader.lines().enumerate() {
            let line = line_result
                .with_context(|| format!("I/O error reading mask BED file: {}", path.display()))?;
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }
            let fields: Vec<&str> = trimmed.split('\t').collect();
            if fields.len() < 3 {
                return Err(anyhow!(
                    "Mask BED line {} has fewer than 3 fields: {}",
                    line_num + 1,
                    trimmed
                ));
            }
            let chr = fields[0].to_string();
            let start_0: u64 = fields[1].parse().with_context(|| {
                format!("Bad start at mask BED line {}: {}", line_num + 1, fields[1])
            })?;
            let end: u64 = fields[2].parse().with_context(|| {
                format!("Bad end at mask BED line {}: {}", line_num + 1, fields[2])
            })?;
            if end == 0 || start_0 >= end {
                continue; // Skip empty intervals silently.
            }
            // Convert to 1-based inclusive.
            let start_1 = start_0 + 1;
            intervals.entry(chr).or_default().push((start_1, end));
        }

        // Sort intervals by start for efficient overlap scanning.
        for v in intervals.values_mut() {
            v.sort_by_key(|&(s, _)| s);
        }

        Ok(MaskIndex { intervals })
    }

    /// Apply masking to a sequence extracted from `chr:region_start-region_end`.
    ///
    /// The sequence string corresponds to 1-based positions
    /// `[region_start, region_end]`. For each overlapping mask interval,
    /// the corresponding bases in the sequence are masked according to `mode`.
    pub fn apply(
        &self,
        seq: &str,
        chr: &str,
        region_start: u64,
        region_end: u64,
        mode: MaskMode,
    ) -> String {
        let Some(chr_intervals) = self.intervals.get(chr) else {
            return seq.to_string();
        };

        let mut bytes = seq.as_bytes().to_vec();

        for &(mask_start, mask_end) in chr_intervals {
            // Skip intervals entirely before or after the region.
            if mask_end < region_start || mask_start > region_end {
                continue;
            }
            // Clamp to the region bounds.
            let overlap_start = mask_start.max(region_start);
            let overlap_end = mask_end.min(region_end);
            // Convert to 0-based index into the sequence string.
            let idx_start = (overlap_start - region_start) as usize;
            let idx_end = (overlap_end - region_start) as usize;

            for b in &mut bytes[idx_start..=idx_end] {
                match mode {
                    MaskMode::Hard => *b = b'N',
                    MaskMode::Soft => *b = b.to_ascii_lowercase(),
                }
            }
        }

        String::from_utf8(bytes).expect("masking produced invalid UTF-8")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn make_index(intervals: &[(&str, u64, u64)]) -> MaskIndex {
        let mut map: HashMap<String, Vec<(u64, u64)>> = HashMap::new();
        for &(chr, start, end) in intervals {
            map.entry(chr.to_string()).or_default().push((start, end));
        }
        for v in map.values_mut() {
            v.sort_by_key(|&(s, _)| s);
        }
        MaskIndex { intervals: map }
    }

    #[test]
    fn hard_mask_full_overlap() {
        let idx = make_index(&[("chr1", 1, 10)]);
        let result = idx.apply("ACGTACGTAC", "chr1", 1, 10, MaskMode::Hard);
        assert_eq!(result, "NNNNNNNNNN");
    }

    #[test]
    fn hard_mask_partial_overlap() {
        let idx = make_index(&[("chr1", 3, 6)]);
        let result = idx.apply("ACGTACGTAC", "chr1", 1, 10, MaskMode::Hard);
        assert_eq!(result, "ACNNNNGTAC");
    }

    #[test]
    fn soft_mask_partial_overlap() {
        let idx = make_index(&[("chr1", 3, 6)]);
        let result = idx.apply("ACGTACGTAC", "chr1", 1, 10, MaskMode::Soft);
        assert_eq!(result, "ACgtacGTAC");
    }

    #[test]
    fn no_overlap() {
        let idx = make_index(&[("chr1", 20, 30)]);
        let result = idx.apply("ACGTACGTAC", "chr1", 1, 10, MaskMode::Hard);
        assert_eq!(result, "ACGTACGTAC");
    }

    #[test]
    fn different_chromosome() {
        let idx = make_index(&[("chr2", 1, 10)]);
        let result = idx.apply("ACGTACGTAC", "chr1", 1, 10, MaskMode::Hard);
        assert_eq!(result, "ACGTACGTAC");
    }

    #[test]
    fn mask_extends_beyond_region() {
        // Mask interval [5, 15] but region is [1, 10]. Should mask [5, 10].
        let idx = make_index(&[("chr1", 5, 15)]);
        let result = idx.apply("ACGTACGTAC", "chr1", 1, 10, MaskMode::Hard);
        assert_eq!(result, "ACGTNNNNNN");
    }

    #[test]
    fn mask_starts_before_region() {
        // Mask interval [1, 5] but region is [3, 12] (10 bases).
        // Overlap is [3, 5], which maps to seq indices 0..=2.
        let idx = make_index(&[("chr1", 1, 5)]);
        let result = idx.apply("ACGTACGTAC", "chr1", 3, 12, MaskMode::Hard);
        assert_eq!(result, "NNNTACGTAC");
    }

    #[test]
    fn multiple_mask_intervals() {
        let idx = make_index(&[("chr1", 1, 3), ("chr1", 8, 10)]);
        let result = idx.apply("ACGTACGTAC", "chr1", 1, 10, MaskMode::Hard);
        assert_eq!(result, "NNNTACGNNN");
    }

    #[test]
    fn from_bed_file() {
        let mut tmp = tempfile::NamedTempFile::new().unwrap();
        // BED: 0-based half-open. chr1:2-5 (0-based) = chr1:3-5 (1-based inclusive).
        writeln!(tmp, "chr1\t2\t5").unwrap();
        writeln!(tmp, "chr1\t7\t9").unwrap();
        tmp.flush().unwrap();
        let idx = MaskIndex::from_bed(tmp.path()).unwrap();
        let result = idx.apply("ACGTACGTAC", "chr1", 1, 10, MaskMode::Hard);
        // Masked positions (1-based): 3-5, 8-9
        assert_eq!(result, "ACNNNCGNNC");
    }

    #[test]
    fn soft_mask_preserves_already_lowercase() {
        let idx = make_index(&[("chr1", 1, 4)]);
        let result = idx.apply("ACgt", "chr1", 1, 4, MaskMode::Soft);
        assert_eq!(result, "acgt");
    }
}
