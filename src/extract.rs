use anyhow::{Result, anyhow};
use std::collections::HashMap;
use std::fs::File;
use std::io::{Read, Seek, SeekFrom};

use crate::fai::FaiRecord;
use crate::region::Region;
#[cfg(test)]
use crate::region::Strand;

/// Extract a region from a FASTA file using the FAI index.
///
/// Uses a single seek and one bulk read for the entire region, then strips
/// newline characters in memory. This minimises syscalls compared to reading
/// line by line.
pub fn extract_region(
    f: &mut File,
    fai: &HashMap<String, FaiRecord>,
    r: &Region,
) -> Result<String> {
    // Coordinates are 1-based; zero is never valid.
    if r.start == 0 {
        return Err(anyhow!(
            "region start must be >= 1 (1-based coordinates), got 0 for '{}'",
            r.chr
        ));
    }

    let rec = fai
        .get(&r.chr)
        .ok_or_else(|| anyhow!("Contig '{}' not in index", r.chr))?;

    if rec.line_bases == 0 {
        return Err(anyhow!(
            "Malformed FAI: line_bases is 0 for contig '{}' (would cause division by zero)",
            r.chr
        ));
    }

    let lb = rec.line_bases;
    let lby = rec.line_bytes;

    // Compute byte range: start offset and end offset (inclusive).
    let start_line = (r.start - 1) / lb;
    let start_col = (r.start - 1) % lb;
    let byte_start = rec.offset + start_line * lby + start_col;

    let end_line = (r.end - 1) / lb;
    let end_col = (r.end - 1) % lb;
    let byte_end = rec.offset + end_line * lby + end_col;

    let read_len = (byte_end - byte_start + 1) as usize;

    // Single seek, single read.
    f.seek(SeekFrom::Start(byte_start))?;
    let mut buf = vec![0u8; read_len];
    f.read_exact(&mut buf)?;

    // Strip newline characters in memory and normalise to uppercase.
    let expected_bases = (r.end - r.start + 1) as usize;
    let mut seq = Vec::<u8>::with_capacity(expected_bases);
    for &b in &buf {
        if b != b'\n' && b != b'\r' {
            seq.push(b.to_ascii_uppercase());
        }
    }

    Ok(String::from_utf8(seq)?)
}

// ─── Unit tests ──────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn extract_region_zero_start_errors() {
        // Build a minimal FASTA file and FAI index for the test.
        let mut tmp = tempfile::NamedTempFile::new().unwrap();
        write!(tmp, ">chr1\nACGTACGT\n").unwrap();
        tmp.flush().unwrap();

        let mut fai = HashMap::new();
        fai.insert(
            "chr1".to_string(),
            FaiRecord {
                length: 8,
                offset: 6,
                line_bases: 8,
                line_bytes: 9,
            },
        );

        let r = Region {
            name: None,
            chr: "chr1".to_string(),
            start: 0,
            end: 5,
            strand: Strand::Unspecified,
        };

        let mut f = File::open(tmp.path()).unwrap();
        let err = extract_region(&mut f, &fai, &r).unwrap_err();
        let msg = err.to_string();
        assert!(
            msg.contains("1-based coordinates"),
            "expected 1-based coordinates message, got: {msg}"
        );
    }
}
