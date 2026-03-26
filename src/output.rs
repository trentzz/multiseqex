use anyhow::{Context, Result, anyhow};
use rayon::prelude::*;
use std::cell::RefCell;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;
use std::sync::Mutex;

use crate::extract::{build_bulk_groups, extract_bulk_group, extract_region, reverse_complement};
use crate::fai::FaiRecord;
use crate::region::Region;

pub fn wrap_fasta(seq: &str, width: usize) -> String {
    if width == 0 {
        return seq.to_string();
    }
    let mut out = String::with_capacity(seq.len() + seq.len() / width + 1);
    for (i, chunk) in seq.as_bytes().chunks(width).enumerate() {
        if i > 0 {
            out.push('\n');
        }
        out.push_str(std::str::from_utf8(chunk).unwrap());
    }
    out
}

fn format_fasta_entry(r: &Region, seq: &str) -> String {
    let header = match &r.name {
        Some(name) => format!(">{name} {}:{}-{}", r.chr, r.start, r.end),
        None => format!(">{}:{}-{}", r.chr, r.start, r.end),
    };
    format!("{header}\n{}\n", wrap_fasta(seq, 60))
}

fn extract_and_format(
    fasta_path: &Path,
    fai_index: &HashMap<String, FaiRecord>,
    r: &Region,
    rc: bool,
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
        let seq = if rc { reverse_complement(&seq) } else { seq };
        let entry = format_fasta_entry(r, &seq);
        Ok((r.clone(), entry))
    })
}

pub fn write_sequences(
    fasta_path: &Path,
    fai_index: &HashMap<String, FaiRecord>,
    regions: &[Region],
    output_file: Option<&Path>,
    output_dir: Option<&Path>,
    is_sv: bool,
    rc: bool,
) -> Result<()> {
    match output_dir {
        Some(dir) => {
            std::fs::create_dir_all(dir)?;
            if is_sv {
                write_sv_per_file_streaming(fasta_path, fai_index, regions, dir, rc)?;
            } else {
                write_per_file_streaming(fasta_path, fai_index, regions, dir, rc)?;
            }
        }
        None => {
            write_streaming_ordered(fasta_path, fai_index, regions, output_file, rc)?;
        }
    }
    Ok(())
}

/// Extract regions in parallel and write them to a single destination
/// (stdout or a file) while preserving the original input order.
///
/// All formatted FASTA entries are buffered in memory before writing, so
/// memory usage scales with total output size rather than with the number
/// of regions alone. For very large extraction jobs where memory is a
/// concern, prefer `--output-dir` which streams each region independently.
fn write_streaming_ordered(
    fasta_path: &Path,
    fai_index: &HashMap<String, FaiRecord>,
    regions: &[Region],
    output_file: Option<&Path>,
    rc: bool,
) -> Result<()> {
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
            let extracted = extract_bulk_group(f, group, regions)?;
            for (orig_idx, seq) in extracted {
                let r = &regions[orig_idx];
                let seq = if rc { reverse_complement(&seq) } else { seq };
                let entry = format_fasta_entry(r, &seq);
                let mut slot = slots[orig_idx].lock().unwrap();
                *slot = Some(entry);
            }
            Ok(())
        })
    })?;
    let mut writer: Box<dyn Write> = match output_file {
        Some(p) => Box::new(BufWriter::new(File::create(p)?)),
        None => Box::new(BufWriter::new(io::stdout())),
    };
    for slot in &slots {
        let guard = slot.lock().unwrap();
        if let Some(entry) = guard.as_ref() {
            writer.write_all(entry.as_bytes())?;
        }
    }
    writer.flush()?;
    Ok(())
}

fn write_per_file_streaming(
    fasta_path: &Path,
    fai_index: &HashMap<String, FaiRecord>,
    regions: &[Region],
    dir: &Path,
    rc: bool,
) -> Result<()> {
    regions.par_iter().try_for_each(|r| -> Result<()> {
        let (region, entry) = extract_and_format(fasta_path, fai_index, r, rc)?;
        let filename = match &region.name {
            Some(name) => format!("{name}_{}_{}.fa", region.start, region.end),
            None => format!("{}_{}_{}.fa", region.chr, region.start, region.end),
        };
        let mut f = File::create(dir.join(filename))?;
        f.write_all(entry.as_bytes())?;
        Ok(())
    })
}

fn write_sv_per_file_streaming(
    fasta_path: &Path,
    fai_index: &HashMap<String, FaiRecord>,
    regions: &[Region],
    dir: &Path,
    rc: bool,
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
            let (r1, entry1) = extract_and_format(fasta_path, fai_index, &pair[0], rc)?;
            let (r2, entry2) = extract_and_format(fasta_path, fai_index, &pair[1], rc)?;
            let filename = match (&r1.name, &r2.name) {
                (Some(n1), Some(n2)) if n1 == n2 => format!(
                    "{n1}_{}_{}_{}_{}_{}_{}.fa",
                    r1.chr, r1.start, r1.end, r2.chr, r2.start, r2.end
                ),
                (Some(n1), Some(n2)) => {
                    return Err(anyhow!("Mismatched names in SV pair: '{n1}' vs '{n2}'"));
                }
                _ => format!(
                    "{}_{}_{}_{}_{}_{}.fa",
                    r1.chr, r1.start, r1.end, r2.chr, r2.start, r2.end
                ),
            };
            let mut f = File::create(dir.join(filename))?;
            f.write_all(entry1.as_bytes())?;
            f.write_all(entry2.as_bytes())?;
            Ok(())
        })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn wrap_fasta_normal() {
        let seq = "A".repeat(120);
        let wrapped = wrap_fasta(&seq, 60);
        let lines: Vec<&str> = wrapped.lines().collect();
        assert_eq!(lines.len(), 2);
        assert_eq!(lines[0].len(), 60);
    }

    #[test]
    fn wrap_fasta_shorter_than_width() {
        assert_eq!(wrap_fasta("ACGT", 60), "ACGT");
    }

    #[test]
    fn wrap_fasta_empty() {
        assert_eq!(wrap_fasta("", 60), "");
    }

    #[test]
    fn wrap_fasta_exact_width() {
        let seq = "A".repeat(60);
        assert_eq!(wrap_fasta(&seq, 60).lines().count(), 1);
    }

    #[test]
    fn wrap_fasta_zero_width() {
        assert_eq!(wrap_fasta("ACGTACGT", 0), "ACGTACGT");
    }

    #[test]
    fn streaming_output_preserves_input_order() {
        let mut tmp = tempfile::NamedTempFile::new().unwrap();
        write!(tmp, ">chr1\nAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCC\n>chr2\nTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAA\n>chr3\nGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTT\n").unwrap();
        tmp.flush().unwrap();
        let mut fai = HashMap::new();
        fai.insert(
            "chr1".to_string(),
            FaiRecord {
                length: 40,
                offset: 6,
                line_bases: 40,
                line_bytes: 41,
            },
        );
        fai.insert(
            "chr2".to_string(),
            FaiRecord {
                length: 40,
                offset: 53,
                line_bases: 40,
                line_bytes: 41,
            },
        );
        fai.insert(
            "chr3".to_string(),
            FaiRecord {
                length: 40,
                offset: 100,
                line_bases: 40,
                line_bytes: 41,
            },
        );
        let regions: Vec<Region> = vec![
            ("chr3", 1u64, 10u64),
            ("chr1", 1, 10),
            ("chr2", 1, 10),
            ("chr1", 11, 20),
            ("chr3", 11, 20),
            ("chr2", 11, 20),
            ("chr1", 21, 30),
            ("chr2", 21, 30),
            ("chr3", 21, 30),
            ("chr1", 31, 40),
            ("chr3", 31, 40),
            ("chr2", 31, 40),
        ]
        .into_iter()
        .map(|(chr, start, end)| Region {
            name: None,
            chr: chr.to_string(),
            start,
            end,
        })
        .collect();
        let out_file = tempfile::NamedTempFile::new().unwrap();
        write_streaming_ordered(tmp.path(), &fai, &regions, Some(out_file.path()), false).unwrap();
        let output = std::fs::read_to_string(out_file.path()).unwrap();
        let headers: Vec<&str> = output.lines().filter(|l| l.starts_with('>')).collect();
        let expected: Vec<String> = regions
            .iter()
            .map(|r| format!(">{}:{}-{}", r.chr, r.start, r.end))
            .collect();
        assert_eq!(headers.len(), expected.len());
        for (i, (got, want)) in headers.iter().zip(expected.iter()).enumerate() {
            assert_eq!(got, want, "order mismatch at {i}");
        }
    }
}
