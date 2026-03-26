use anyhow::{Context, Result, anyhow};
use rayon::prelude::*;
use std::cell::RefCell;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

use crate::extract::extract_region;
use crate::fai::FaiRecord;
use crate::region::Region;

/// Wrap a sequence string to a fixed line width for FASTA output.
/// If `width` is 0, the sequence is returned unwrapped.
pub fn wrap_fasta(seq: &str, width: usize) -> String {
    if width == 0 {
        return seq.to_string();
    }
    let mut out = String::with_capacity(seq.len() + seq.len() / width + 1);
    for (i, chunk) in seq.as_bytes().chunks(width).enumerate() {
        if i > 0 {
            out.push('\n');
        }
        // SAFETY: input is valid UTF-8 ASCII, so chunks are too.
        out.push_str(std::str::from_utf8(chunk).unwrap());
    }
    out
}

/// Format a single region's extraction as a FASTA entry string.
fn format_fasta_entry(r: &Region, seq: &str) -> String {
    let header = match &r.name {
        Some(name) => format!(">{name} {}:{}-{}", r.chr, r.start, r.end),
        None => format!(">{}:{}-{}", r.chr, r.start, r.end),
    };
    format!("{header}\n{}\n", wrap_fasta(seq, 60))
}

/// Extract a single region and return the formatted FASTA entry.
///
/// Uses a thread-local file handle for reuse across calls on the same thread.
fn extract_and_format(
    fasta_path: &Path,
    fai_index: &HashMap<String, FaiRecord>,
    r: &Region,
) -> Result<(Region, String)> {
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
        let seq = extract_region(f, fai_index, r)?;
        let entry = format_fasta_entry(r, &seq);
        Ok((r.clone(), entry))
    })
}

/// Extract and write all sequences to the requested destination.
///
/// For `--output-dir` mode, files are written as each extraction completes,
/// avoiding the need to hold all sequences in memory simultaneously.
/// For single-file and stdout output, results are collected to preserve
/// deterministic (input) order.
pub fn write_sequences(
    fasta_path: &Path,
    fai_index: &HashMap<String, FaiRecord>,
    regions: &[Region],
    output_file: Option<&Path>,
    output_dir: Option<&Path>,
    is_sv: bool,
) -> Result<()> {
    match output_dir {
        Some(dir) => {
            std::fs::create_dir_all(dir)?;
            if is_sv {
                write_sv_per_file_streaming(fasta_path, fai_index, regions, dir)?;
            } else {
                write_per_file_streaming(fasta_path, fai_index, regions, dir)?;
            }
        }
        None => {
            // Collect all results to maintain deterministic input order.
            let sequences: Vec<(Region, String)> = regions
                .par_iter()
                .map(|r| extract_and_format(fasta_path, fai_index, r))
                .collect::<Result<_>>()?;

            let mut writer: Box<dyn Write> = match output_file {
                Some(p) => Box::new(BufWriter::new(File::create(p)?)),
                None => Box::new(BufWriter::new(io::stdout())),
            };
            for (_, entry) in &sequences {
                writer.write_all(entry.as_bytes())?;
            }
            writer.flush()?;
        }
    }
    Ok(())
}

/// Stream one FASTA file per region directly to disk, without collecting all
/// sequences into memory first.
fn write_per_file_streaming(
    fasta_path: &Path,
    fai_index: &HashMap<String, FaiRecord>,
    regions: &[Region],
    dir: &Path,
) -> Result<()> {
    regions.par_iter().try_for_each(|r| -> Result<()> {
        let (region, entry) = extract_and_format(fasta_path, fai_index, r)?;
        let filename = match &region.name {
            Some(name) => format!("{name}_{}_{}.fa", region.start, region.end),
            None => format!("{}_{}_{}.fa", region.chr, region.start, region.end),
        };
        let mut f = File::create(dir.join(filename))?;
        f.write_all(entry.as_bytes())?;
        Ok(())
    })
}

/// Stream one FASTA file per SV pair directly to disk.
fn write_sv_per_file_streaming(
    fasta_path: &Path,
    fai_index: &HashMap<String, FaiRecord>,
    regions: &[Region],
    dir: &Path,
) -> Result<()> {
    if !regions.len().is_multiple_of(2) {
        return Err(anyhow!(
            "SV extraction produced {} regions (expected even number for pairs)",
            regions.len()
        ));
    }

    regions
        .chunks(2)
        .collect::<Vec<_>>()
        .par_iter()
        .try_for_each(|pair| -> Result<()> {
            let (r1, entry1) = extract_and_format(fasta_path, fai_index, &pair[0])?;
            let (r2, entry2) = extract_and_format(fasta_path, fai_index, &pair[1])?;

            let filename = match (&r1.name, &r2.name) {
                (Some(n1), Some(n2)) if n1 == n2 => {
                    format!(
                        "{n1}_{}_{}_{}_{}_{}_{}.fa",
                        r1.chr, r1.start, r1.end, r2.chr, r2.start, r2.end
                    )
                }
                (Some(n1), Some(n2)) => {
                    return Err(anyhow!("Mismatched names in SV pair: '{n1}' vs '{n2}'"));
                }
                _ => {
                    format!(
                        "{}_{}_{}_{}_{}_{}.fa",
                        r1.chr, r1.start, r1.end, r2.chr, r2.start, r2.end
                    )
                }
            };

            let mut f = File::create(dir.join(filename))?;
            f.write_all(entry1.as_bytes())?;
            f.write_all(entry2.as_bytes())?;
            Ok(())
        })
}

// ─── Unit tests ──────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── wrap_fasta ───────────────────────────────────────────────────────

    #[test]
    fn wrap_fasta_normal() {
        let seq = "A".repeat(120);
        let wrapped = wrap_fasta(&seq, 60);
        let lines: Vec<&str> = wrapped.lines().collect();
        assert_eq!(lines.len(), 2);
        assert_eq!(lines[0].len(), 60);
        assert_eq!(lines[1].len(), 60);
    }

    #[test]
    fn wrap_fasta_shorter_than_width() {
        let wrapped = wrap_fasta("ACGT", 60);
        assert_eq!(wrapped, "ACGT");
    }

    #[test]
    fn wrap_fasta_empty() {
        let wrapped = wrap_fasta("", 60);
        assert_eq!(wrapped, "");
    }

    #[test]
    fn wrap_fasta_exact_width() {
        let seq = "A".repeat(60);
        let wrapped = wrap_fasta(&seq, 60);
        assert_eq!(wrapped.lines().count(), 1);
    }

    #[test]
    fn wrap_fasta_zero_width() {
        let wrapped = wrap_fasta("ACGTACGT", 0);
        assert_eq!(wrapped, "ACGTACGT");
    }
}
