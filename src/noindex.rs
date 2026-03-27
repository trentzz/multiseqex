//! Sequential scan mode for FASTA extraction without an FAI index.
//!
//! When `--no-index` is set, this module reads the FASTA sequentially,
//! building an in-memory representation, and extracts the requested regions.
//! This supports piped input (stdin via `-` as the FASTA path).
//!
//! For truly large files, this loads the entire FASTA into memory. For those
//! cases, a pre-built FAI index is recommended.

use anyhow::{Result, anyhow};
use std::collections::HashMap;
use std::io::{BufRead, BufReader, Read};

use crate::fai::FaiRecord;
use crate::region::Region;

/// An in-memory contig: name and sequence bytes (uppercase, no newlines).
struct InMemoryContig {
    sequence: Vec<u8>,
}

/// The result of loading a FASTA into memory: sequence data and FAI records.
pub type InMemoryFasta = (HashMap<String, Vec<u8>>, HashMap<String, FaiRecord>);

/// Read a FASTA from any reader and build an in-memory index.
///
/// Returns a map of contig name to sequence, plus a FaiRecord map for
/// compatibility with the validation and extraction pipeline.
pub fn load_fasta_into_memory(reader: Box<dyn Read>) -> Result<InMemoryFasta> {
    let buf_reader = BufReader::new(reader);
    let mut contigs: HashMap<String, InMemoryContig> = HashMap::new();
    let mut fai: HashMap<String, FaiRecord> = HashMap::new();
    let mut current_name: Option<String> = None;
    let mut current_seq: Vec<u8> = Vec::new();

    for line_result in buf_reader.lines() {
        let line = line_result?;
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }

        if let Some(header) = trimmed.strip_prefix('>') {
            // Flush previous contig.
            if let Some(name) = current_name.take() {
                let length = current_seq.len() as u64;
                fai.insert(
                    name.clone(),
                    FaiRecord {
                        length,
                        offset: 0,     // Not meaningful for in-memory mode.
                        line_bases: 0, // Not meaningful for in-memory mode.
                        line_bytes: 0, // Not meaningful for in-memory mode.
                    },
                );
                contigs.insert(
                    name,
                    InMemoryContig {
                        sequence: std::mem::take(&mut current_seq),
                    },
                );
            }

            let name = header.split_whitespace().next().unwrap_or("").to_string();
            if name.is_empty() {
                return Err(anyhow!("FASTA header with empty contig name"));
            }
            current_name = Some(name);
            current_seq.clear();
        } else {
            // Sequence line: accumulate uppercase bases.
            for &b in trimmed.as_bytes() {
                if b.is_ascii_alphabetic() {
                    current_seq.push(b.to_ascii_uppercase());
                }
            }
        }
    }

    // Flush the last contig.
    if let Some(name) = current_name.take() {
        let length = current_seq.len() as u64;
        fai.insert(
            name.clone(),
            FaiRecord {
                length,
                offset: 0,
                line_bases: 0,
                line_bytes: 0,
            },
        );
        contigs.insert(
            name,
            InMemoryContig {
                sequence: current_seq,
            },
        );
    }

    let sequences: HashMap<String, Vec<u8>> = contigs
        .into_iter()
        .map(|(name, c)| (name, c.sequence))
        .collect();

    Ok((sequences, fai))
}

/// Extract a region from the in-memory sequence store.
pub fn extract_region_from_memory(
    sequences: &HashMap<String, Vec<u8>>,
    r: &Region,
) -> Result<String> {
    let seq = sequences
        .get(&r.chr)
        .ok_or_else(|| anyhow!("Contig '{}' not found in FASTA", r.chr))?;

    if r.start == 0 {
        return Err(anyhow!(
            "Region start must be >= 1 (1-based coordinates), got 0 for '{}'",
            r.chr
        ));
    }

    let start_idx = (r.start - 1) as usize;
    let end_idx = (r.end as usize).min(seq.len());

    if start_idx >= seq.len() {
        return Err(anyhow!(
            "Region {}:{}-{} is beyond the contig length ({})",
            r.chr,
            r.start,
            r.end,
            seq.len()
        ));
    }

    let slice = &seq[start_idx..end_idx];
    Ok(String::from_utf8(slice.to_vec())?)
}

// ─── Unit tests ──────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn load_simple_fasta() {
        let fasta = b">chr1\nACGTACGT\n>chr2\nTTTTGGGG\n";
        let reader: Box<dyn Read> = Box::new(Cursor::new(fasta.to_vec()));
        let (seqs, fai) = load_fasta_into_memory(reader).unwrap();

        assert_eq!(seqs.len(), 2);
        assert_eq!(seqs["chr1"], b"ACGTACGT");
        assert_eq!(seqs["chr2"], b"TTTTGGGG");
        assert_eq!(fai["chr1"].length, 8);
        assert_eq!(fai["chr2"].length, 8);
    }

    #[test]
    fn load_multiline_fasta() {
        let fasta = b">chr1\nACGT\nACGT\n";
        let reader: Box<dyn Read> = Box::new(Cursor::new(fasta.to_vec()));
        let (seqs, _) = load_fasta_into_memory(reader).unwrap();

        assert_eq!(seqs["chr1"], b"ACGTACGT");
    }

    #[test]
    fn extract_from_memory() {
        let fasta = b">chr1\nACGTACGTACGT\n";
        let reader: Box<dyn Read> = Box::new(Cursor::new(fasta.to_vec()));
        let (seqs, _) = load_fasta_into_memory(reader).unwrap();

        let r = Region {
            name: None,
            chr: "chr1".to_string(),
            start: 3,
            end: 8,
            strand: None,
        };
        let seq = extract_region_from_memory(&seqs, &r).unwrap();
        assert_eq!(seq, "GTACGT");
    }

    #[test]
    fn extract_clamped_to_length() {
        let fasta = b">chr1\nACGT\n";
        let reader: Box<dyn Read> = Box::new(Cursor::new(fasta.to_vec()));
        let (seqs, _) = load_fasta_into_memory(reader).unwrap();

        let r = Region {
            name: None,
            chr: "chr1".to_string(),
            start: 2,
            end: 100,
            strand: None,
        };
        let seq = extract_region_from_memory(&seqs, &r).unwrap();
        assert_eq!(seq, "CGT");
    }

    #[test]
    fn extract_unknown_contig_errors() {
        let fasta = b">chr1\nACGT\n";
        let reader: Box<dyn Read> = Box::new(Cursor::new(fasta.to_vec()));
        let (seqs, _) = load_fasta_into_memory(reader).unwrap();

        let r = Region {
            name: None,
            chr: "chrZ".to_string(),
            start: 1,
            end: 4,
            strand: None,
        };
        assert!(extract_region_from_memory(&seqs, &r).is_err());
    }

    #[test]
    fn empty_fasta() {
        let fasta = b"";
        let reader: Box<dyn Read> = Box::new(Cursor::new(fasta.to_vec()));
        let (seqs, fai) = load_fasta_into_memory(reader).unwrap();
        assert!(seqs.is_empty());
        assert!(fai.is_empty());
    }
}
