use anyhow::{Result, anyhow};
use std::cmp::min;
use std::collections::HashMap;
use std::fs::File;
use std::io::{Read, Seek, SeekFrom};

use crate::fai::FaiRecord;
use crate::region::Region;

/// Extract a region from a FASTA file using the FAI index.
/// Accepts a mutable file handle to allow reuse across calls on the same thread.
pub(crate) fn extract_region(
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

    let mut seq = Vec::<u8>::with_capacity((r.end - r.start + 1) as usize);
    let mut pos = r.start;

    while pos <= r.end {
        let line_idx = (pos - 1) / lb;
        let in_line_offset = (pos - 1) % lb;
        let run = min(lb - in_line_offset, r.end - pos + 1);
        let byte_pos = rec.offset + line_idx * lby + in_line_offset;

        f.seek(SeekFrom::Start(byte_pos))?;
        let mut buf = vec![0u8; run as usize];
        f.read_exact(&mut buf)?;
        seq.extend_from_slice(&buf);

        pos += run;
    }

    // Normalise to uppercase.
    for b in &mut seq {
        b.make_ascii_uppercase();
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
