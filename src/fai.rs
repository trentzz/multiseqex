use anyhow::{Context, Result};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

/// One record from a `.fai` index file.
#[derive(Debug, Clone)]
pub(crate) struct FaiRecord {
    /// Total number of bases in this contig.
    pub(crate) length: u64,
    /// Byte offset of the first base in the FASTA file.
    pub(crate) offset: u64,
    /// Number of sequence bases per line.
    pub(crate) line_bases: u64,
    /// Number of bytes per line (bases + newline characters).
    pub(crate) line_bytes: u64,
}

/// Compute the `.fai` path for a given FASTA file.
pub(crate) fn fai_path_for(fasta: &Path) -> PathBuf {
    let mut s = fasta.as_os_str().to_owned();
    s.push(".fai");
    PathBuf::from(s)
}

/// Warn if the FAI file is older than the FASTA file, suggesting it may be stale.
pub(crate) fn check_fai_staleness(fasta: &Path, fai: &Path) {
    let fasta_mtime = std::fs::metadata(fasta).and_then(|m| m.modified());
    let fai_mtime = std::fs::metadata(fai).and_then(|m| m.modified());

    if let (Ok(fasta_t), Ok(fai_t)) = (fasta_mtime, fai_mtime)
        && fai_t < fasta_t
    {
        eprintln!(
            "Warning: FAI index '{}' is older than FASTA '{}'. \
             The index may be stale. Consider rebuilding it.",
            fai.display(),
            fasta.display()
        );
    }
}

/// Build a minimal `.fai` index from a FASTA file.
pub(crate) fn build_fai(fasta: &Path, fai_out: &Path) -> Result<()> {
    let f = File::open(fasta)
        .with_context(|| format!("Cannot open FASTA for indexing: {}", fasta.display()))?;
    let mut reader = BufReader::new(f);
    let mut out = BufWriter::new(
        File::create(fai_out)
            .with_context(|| format!("Cannot create FAI: {}", fai_out.display()))?,
    );

    let mut pos: u64 = 0;
    let mut line = String::new();

    let mut current_name: Option<String> = None;
    let mut seq_offset: u64 = 0;
    let mut seq_len: u64 = 0;
    let mut line_bases: u64 = 0;
    let mut line_bytes: u64 = 0;
    let mut first_seq_line = true;

    loop {
        line.clear();
        let n = reader.read_line(&mut line)?;
        if n == 0 {
            // EOF: flush last contig.
            if let Some(name) = current_name.take() {
                writeln!(
                    out,
                    "{name}\t{seq_len}\t{seq_offset}\t{line_bases}\t{line_bytes}"
                )?;
            }
            break;
        }

        let raw = line.as_bytes();
        let linelen = raw.len() as u64;

        if raw.starts_with(b">") {
            // New contig header: flush the previous one.
            let header_name = parse_fasta_header(&line);
            if header_name.is_empty() {
                eprintln!(
                    "Warning: empty contig name at byte offset {}. Header line is bare '>' or '>  '.",
                    pos
                );
            }
            if let Some(name) = current_name.replace(header_name) {
                writeln!(
                    out,
                    "{name}\t{seq_len}\t{seq_offset}\t{line_bases}\t{line_bytes}"
                )?;
            }
            seq_offset = pos + linelen;
            seq_len = 0;
            line_bases = 0;
            line_bytes = 0;
            first_seq_line = true;
        } else {
            let bases = count_bases(raw);
            seq_len += bases;
            if first_seq_line {
                line_bases = bases;
                line_bytes = linelen;
                first_seq_line = false;
            } else if bases > 0 && bases != line_bases && linelen != line_bytes {
                // Non-final lines with a different width indicate an inconsistent FASTA.
                // Only warn if this is not the last (possibly shorter) line of the contig.
                // We cannot know for certain it is non-final here, so we check that
                // the line is shorter than expected (final lines are allowed to differ).
                if linelen >= line_bytes {
                    eprintln!(
                        "Warning: inconsistent line width in contig '{}': expected {} bases/{} bytes, got {} bases/{} bytes",
                        current_name.as_deref().unwrap_or("?"),
                        line_bases,
                        line_bytes,
                        bases,
                        linelen
                    );
                }
            }
        }
        pos += linelen;
    }

    out.flush()?;
    Ok(())
}

/// Extract the first whitespace-delimited token after `>`.
pub(crate) fn parse_fasta_header(s: &str) -> String {
    s.trim_start_matches('>')
        .split_whitespace()
        .next()
        .unwrap_or("")
        .to_string()
}

/// Count ASCII alphabetic characters in a raw line.
pub(crate) fn count_bases(raw: &[u8]) -> u64 {
    raw.iter().filter(|b| b.is_ascii_alphabetic()).count() as u64
}

/// Read a `.fai` file into a contig-name to `FaiRecord` map.
pub(crate) fn read_fai(fai_path: &Path) -> Result<HashMap<String, FaiRecord>> {
    let f =
        File::open(fai_path).with_context(|| format!("Cannot open FAI: {}", fai_path.display()))?;
    let reader = BufReader::new(f);
    let mut map = HashMap::new();

    for (i, line) in reader.lines().enumerate() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 5 {
            return Err(anyhow::anyhow!("Malformed FAI line {}: {}", i + 1, line));
        }
        let contig_name = parts[0].to_string();
        if map.contains_key(&contig_name) {
            eprintln!(
                "Warning: duplicate contig name '{}' in FAI (line {}). Earlier entry will be overwritten.",
                contig_name,
                i + 1
            );
        }
        map.insert(
            contig_name,
            FaiRecord {
                length: parts[1]
                    .parse()
                    .with_context(|| format!("Bad length, FAI line {}", i + 1))?,
                offset: parts[2]
                    .parse()
                    .with_context(|| format!("Bad offset, FAI line {}", i + 1))?,
                line_bases: parts[3]
                    .parse()
                    .with_context(|| format!("Bad line_bases, FAI line {}", i + 1))?,
                line_bytes: parts[4]
                    .parse()
                    .with_context(|| format!("Bad line_bytes, FAI line {}", i + 1))?,
            },
        );
    }
    Ok(map)
}

// ─── Unit tests ──────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── count_bases ──────────────────────────────────────────────────────

    #[test]
    fn count_bases_normal() {
        assert_eq!(count_bases(b"ACGTNN\n"), 6);
    }

    #[test]
    fn count_bases_empty() {
        assert_eq!(count_bases(b""), 0);
    }

    #[test]
    fn count_bases_non_alpha() {
        assert_eq!(count_bases(b"123\n"), 0);
    }

    #[test]
    fn count_bases_mixed() {
        assert_eq!(count_bases(b"AC1GT\r\n"), 4);
    }

    // ── parse_fasta_header ───────────────────────────────────────────────

    #[test]
    fn parse_fasta_header_normal() {
        assert_eq!(parse_fasta_header(">chr1"), "chr1");
    }

    #[test]
    fn parse_fasta_header_with_description() {
        assert_eq!(parse_fasta_header(">chr1 some description"), "chr1");
    }

    #[test]
    fn parse_fasta_header_empty() {
        assert_eq!(parse_fasta_header(">"), "");
    }

    #[test]
    fn parse_fasta_header_whitespace_only() {
        assert_eq!(parse_fasta_header(">  "), "");
    }
}
