use anyhow::{Context, Result, anyhow};
use std::collections::HashMap;
use std::fs::File;
use std::io::Read;
use std::path::{Path, PathBuf};

use crate::fai::FaiRecord;
use crate::region::Region;

/// Check whether a file starts with gzip magic bytes (0x1f 0x8b).
pub fn is_gzip(fasta: &Path) -> Result<bool> {
    let mut f =
        File::open(fasta).with_context(|| format!("Cannot open FASTA: {}", fasta.display()))?;
    let mut magic = [0u8; 2];
    if f.read_exact(&mut magic).is_ok() && magic == [0x1f, 0x8b] {
        return Ok(true);
    }
    Ok(false)
}

/// Decompress a gzip/bgzip FASTA to a temporary file.
///
/// Returns the path to the decompressed temp file. The caller must keep the
/// returned `TempDir` alive for the duration of use, since dropping it deletes
/// the temporary directory and its contents.
///
/// This is a pragmatic approach: for large bgzipped FASTAs, pre-decompression
/// with `bgzip -d` or `gunzip` is recommended.
pub fn decompress_gzip_to_temp(fasta: &Path) -> Result<(PathBuf, tempfile::TempDir)> {
    let f = File::open(fasta).with_context(|| format!("Cannot open FASTA: {}", fasta.display()))?;
    let mut decoder = flate2::read::GzDecoder::new(f);

    let tmp_dir = tempfile::TempDir::new()
        .with_context(|| "Failed to create temporary directory for bgzip decompression")?;
    let decompressed_path = tmp_dir.path().join("decompressed.fa");
    let mut out = File::create(&decompressed_path)
        .with_context(|| format!("Cannot create temp file: {}", decompressed_path.display()))?;

    std::io::copy(&mut decoder, &mut out)
        .with_context(|| format!("Failed to decompress: {}", fasta.display()))?;

    Ok((decompressed_path, tmp_dir))
}

/// Check whether a FASTA file is gzip/bgzip compressed. If so, decompress it
/// to a temporary file and return the decompressed path. If not, return the
/// original path unchanged.
///
/// When decompression occurs, the returned `Option<TempDir>` holds the temp
/// directory. The caller must keep it alive until extraction is complete.
pub fn resolve_bgzip(fasta: &Path, quiet: bool) -> Result<(PathBuf, Option<tempfile::TempDir>)> {
    if !is_gzip(fasta)? {
        return Ok((fasta.to_path_buf(), None));
    }

    if !quiet {
        eprintln!(
            "Detected gzip/bgzip compressed FASTA: {}. Decompressing to a temporary file.",
            fasta.display()
        );
        eprintln!(
            "Note: for large files, pre-decompression with 'bgzip -d' or 'gunzip' is recommended."
        );
    }

    let (decompressed, tmp_dir) = decompress_gzip_to_temp(fasta)?;
    Ok((decompressed, Some(tmp_dir)))
}

/// Reject gzip/bgzip compressed files by checking magic bytes (0x1f 0x8b).
pub fn detect_gzip_and_reject(fasta: &Path) -> Result<()> {
    if is_gzip(fasta)? {
        return Err(anyhow!(
            "File appears to be gzip/bgzip compressed: {}. Decompress it first (e.g. gunzip or bgzip -d).",
            fasta.display()
        ));
    }
    Ok(())
}

/// Validate that all regions reference known contigs and clamp to contig bounds.
pub fn validate_and_clamp_regions(
    regions: &mut [Region],
    fai: &HashMap<String, FaiRecord>,
) -> Result<()> {
    let mut missing = Vec::new();

    for r in regions.iter_mut() {
        if let Some(rec) = fai.get(&r.chr) {
            r.start = r.start.max(1).min(rec.length);
            r.end = r.end.max(1).min(rec.length);
            if r.start > r.end {
                std::mem::swap(&mut r.start, &mut r.end);
            }
        } else {
            missing.push(r.chr.clone());
        }
    }

    if !missing.is_empty() {
        missing.sort();
        missing.dedup();
        return Err(anyhow!(
            "Contigs not found in FASTA/FAI: {}",
            missing.join(", ")
        ));
    }
    Ok(())
}

// ─── Unit tests ──────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── detect_gzip_and_reject ───────────────────────────────────────────

    #[test]
    fn detect_gzip_rejects_gzip_file() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(tmp.path(), [0x1f, 0x8b, 0x08, 0x00]).unwrap();
        assert!(detect_gzip_and_reject(tmp.path()).is_err());
    }

    #[test]
    fn detect_gzip_accepts_plain_fasta() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(tmp.path(), b">chr1\nACGT\n").unwrap();
        assert!(detect_gzip_and_reject(tmp.path()).is_ok());
    }

    #[test]
    fn detect_gzip_accepts_empty_file() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        // Empty file: read_exact fails, so no rejection.
        assert!(detect_gzip_and_reject(tmp.path()).is_ok());
    }
}
