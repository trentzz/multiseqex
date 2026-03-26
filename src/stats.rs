//! Statistics mode: compute per-region summary statistics instead of
//! extracting sequences.
//!
//! Outputs a TSV table with columns: chr, start, end, name, length,
//! gc_percent, n_count, masked_count.

use anyhow::{Context, Result};
use rayon::prelude::*;
use std::cell::RefCell;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::sync::Mutex;

use crate::extract::build_bulk_groups;
use crate::fai::FaiRecord;
use crate::region::Region;

/// Statistics for a single region.
pub struct RegionStats {
    pub length: u64,
    pub gc_count: u64,
    pub n_count: u64,
    pub masked_count: u64,
}

impl RegionStats {
    /// GC percentage as a float (0.0 to 100.0).
    pub fn gc_percent(&self) -> f64 {
        if self.length == 0 {
            return 0.0;
        }
        (self.gc_count as f64 / self.length as f64) * 100.0
    }
}

/// Compute statistics for a sequence string.
///
/// The sequence is expected to contain only DNA bases (A, C, G, T, N) in
/// upper or lower case. Lowercase bases are counted as masked (soft-masked).
pub fn compute_stats(seq: &str) -> RegionStats {
    let mut gc_count: u64 = 0;
    let mut n_count: u64 = 0;
    let mut masked_count: u64 = 0;
    let length = seq.len() as u64;

    for &b in seq.as_bytes() {
        match b {
            b'G' | b'C' => gc_count += 1,
            b'g' | b'c' => {
                gc_count += 1;
                masked_count += 1;
            }
            b'N' => n_count += 1,
            b'n' => {
                n_count += 1;
                masked_count += 1;
            }
            b'a' | b't' => masked_count += 1,
            _ => {}
        }
    }

    RegionStats {
        length,
        gc_count,
        n_count,
        masked_count,
    }
}

/// Extract sequences without uppercasing, preserving case for masked-base detection.
///
/// The standard `extract_region` uppercases everything. For stats mode we need
/// the original case, so we re-extract from raw bytes here.
fn extract_raw_sequence(
    f: &mut File,
    fai: &HashMap<String, FaiRecord>,
    r: &Region,
) -> Result<String> {
    use std::io::{Read, Seek, SeekFrom};

    let rec = fai
        .get(&r.chr)
        .ok_or_else(|| anyhow::anyhow!("Contig '{}' not in index", r.chr))?;
    let lb = rec.line_bases;
    let lby = rec.line_bytes;
    let start_line = (r.start - 1) / lb;
    let start_col = (r.start - 1) % lb;
    let byte_start = rec.offset + start_line * lby + start_col;
    let end_line = (r.end - 1) / lb;
    let end_col = (r.end - 1) % lb;
    let byte_end = rec.offset + end_line * lby + end_col;
    let read_len = (byte_end - byte_start + 1) as usize;
    f.seek(SeekFrom::Start(byte_start))?;
    let mut buf = vec![0u8; read_len];
    f.read_exact(&mut buf)?;
    // Strip newlines but preserve case.
    let seq: Vec<u8> = buf
        .into_iter()
        .filter(|&b| b != b'\n' && b != b'\r')
        .collect();
    Ok(String::from_utf8(seq)?)
}

/// Write statistics for all regions as a TSV table.
pub fn write_stats(
    fasta_path: &std::path::Path,
    fai_index: &HashMap<String, FaiRecord>,
    regions: &[Region],
) -> Result<()> {
    // We need raw (case-preserving) extraction for masked-base counting.
    // Use parallel extraction via bulk groups, but with raw sequences.
    let groups = build_bulk_groups(regions, fai_index)?;
    let slots: Vec<Mutex<Option<String>>> = (0..regions.len()).map(|_| Mutex::new(None)).collect();

    groups.par_iter().try_for_each(|group| -> Result<()> {
        thread_local! {
            static TL_FILE: RefCell<Option<File>> = const { RefCell::new(None) };
        }
        TL_FILE.with(|cell| {
            let mut borrow = cell.borrow_mut();
            if borrow.is_none() {
                *borrow = Some(
                    File::open(fasta_path)
                        .with_context(|| format!("Cannot open FASTA: {}", fasta_path.display()))?,
                );
            }
            let f = borrow.as_mut().unwrap();
            // Extract raw sequences for this group.
            // We reuse bulk group structure but do raw extraction per region.
            for &orig_idx in &group.indices {
                let r = &regions[orig_idx];
                let seq = extract_raw_sequence(f, fai_index, r)?;
                let mut slot = slots[orig_idx].lock().unwrap();
                *slot = Some(seq);
            }
            Ok(())
        })
    })?;

    let mut writer = BufWriter::new(io::stdout());
    // Header row.
    writeln!(
        writer,
        "chr\tstart\tend\tname\tlength\tgc_percent\tn_count\tmasked_count"
    )?;

    for (i, r) in regions.iter().enumerate() {
        let guard = slots[i].lock().unwrap();
        let seq = guard.as_ref().expect("sequence not extracted");
        let stats = compute_stats(seq);
        let name = r.name.as_deref().unwrap_or(".");
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{:.2}\t{}\t{}",
            r.chr,
            r.start,
            r.end,
            name,
            stats.length,
            stats.gc_percent(),
            stats.n_count,
            stats.masked_count
        )?;
    }
    writer.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stats_all_gc() {
        let s = compute_stats("GCGCGC");
        assert_eq!(s.length, 6);
        assert_eq!(s.gc_count, 6);
        assert_eq!(s.n_count, 0);
        assert_eq!(s.masked_count, 0);
        assert!((s.gc_percent() - 100.0).abs() < 0.01);
    }

    #[test]
    fn stats_mixed() {
        let s = compute_stats("ACGT");
        assert_eq!(s.length, 4);
        assert_eq!(s.gc_count, 2);
        assert_eq!(s.n_count, 0);
        assert_eq!(s.masked_count, 0);
        assert!((s.gc_percent() - 50.0).abs() < 0.01);
    }

    #[test]
    fn stats_with_n() {
        let s = compute_stats("ACNGT");
        assert_eq!(s.length, 5);
        assert_eq!(s.gc_count, 2);
        assert_eq!(s.n_count, 1);
        assert_eq!(s.masked_count, 0);
    }

    #[test]
    fn stats_masked() {
        let s = compute_stats("acgt");
        assert_eq!(s.length, 4);
        assert_eq!(s.gc_count, 2);
        assert_eq!(s.n_count, 0);
        assert_eq!(s.masked_count, 4);
    }

    #[test]
    fn stats_mixed_case() {
        // A=upper, C=upper gc, g=lower gc+masked, t=lower masked,
        // N=upper n, n=lower n+masked.
        let s = compute_stats("ACgtNn");
        assert_eq!(s.length, 6);
        assert_eq!(s.gc_count, 2); // C, g
        assert_eq!(s.n_count, 2); // N, n
        assert_eq!(s.masked_count, 3); // g, t, n
    }

    #[test]
    fn stats_empty() {
        let s = compute_stats("");
        assert_eq!(s.length, 0);
        assert!((s.gc_percent() - 0.0).abs() < 0.01);
    }
}
