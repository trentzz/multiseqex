use anyhow::{Context, Result, anyhow};
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use std::cell::RefCell;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::Mutex;

use crate::extract::{build_bulk_groups, extract_bulk_group, extract_region, reverse_complement};
use crate::fai::FaiRecord;
use crate::mask::{MaskIndex, MaskMode};
use crate::region::Region;
use crate::template::expand_template;
use crate::transform::TransformConfig;

/// Output format configuration.
#[derive(Debug, Clone, Default)]
pub struct OutputConfig {
    /// Emit FASTQ instead of FASTA.
    pub fastq: bool,
    /// Quality character for FASTQ output (default 'I', phred 40).
    pub qual_char: char,
    /// Optional name template for headers.
    pub name_template: Option<String>,
    /// Optional VCF description to append to headers (per region, indexed by position).
    pub vcf_descriptions: Vec<String>,
    /// Optional mask index for masking extracted sequences.
    pub mask_index: Option<MaskIndex>,
    /// Masking mode (hard or soft). Only used when mask_index is Some.
    pub mask_mode: MaskMode,
    /// Sequence transform configuration.
    pub transform: TransformConfig,
}

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

/// Format the strand suffix for a FASTA header.
fn strand_suffix(r: &Region) -> &'static str {
    match r.strand {
        Some('+') => "(+)",
        Some('-') => "(-)",
        Some('.') => "(.)",
        _ => "",
    }
}

#[cfg(test)]
fn format_fasta_entry(r: &Region, seq: &str, line_width: usize) -> String {
    let suffix = strand_suffix(r);
    let header = match &r.name {
        Some(name) => format!(">{name} {}:{}-{}{suffix}", r.chr, r.start, r.end),
        None => format!(">{}:{}-{}{suffix}", r.chr, r.start, r.end),
    };
    format!("{header}\n{}\n", wrap_fasta(seq, line_width))
}

/// Format a FASTA/FASTQ entry with optional template and description.
fn format_entry_with_config(
    r: &Region,
    seq: &str,
    line_width: usize,
    index: usize,
    config: &OutputConfig,
) -> String {
    let description = if index < config.vcf_descriptions.len() {
        Some(config.vcf_descriptions[index].as_str())
    } else {
        None
    };

    let header_content = if let Some(template) = &config.name_template {
        expand_template(template, r, index + 1) // 1-based index
    } else {
        let suffix = strand_suffix(r);
        match &r.name {
            Some(name) => format!("{name} {}:{}-{}{suffix}", r.chr, r.start, r.end),
            None => format!("{}:{}-{}{suffix}", r.chr, r.start, r.end),
        }
    };

    let full_header = match description {
        Some(desc) => format!("{header_content} {desc}"),
        None => header_content,
    };

    if config.fastq {
        let qual: String = std::iter::repeat_n(config.qual_char, seq.len()).collect();
        format!("@{full_header}\n{seq}\n+\n{qual}\n")
    } else {
        format!(">{full_header}\n{}\n", wrap_fasta(seq, line_width))
    }
}

/// Format a region as a TSV line: chr, start, end, name, sequence.
fn format_tab_entry(r: &Region, seq: &str) -> String {
    let name = r.name.as_deref().unwrap_or(".");
    format!("{}\t{}\t{}\t{}\t{}\n", r.chr, r.start, r.end, name, seq)
}

/// Apply per-region strand reverse complement and global --rc flag.
///
/// When a region has strand == Some('-'), we reverse-complement it first.
/// If --rc is also set, the two cancel out (- strand + --rc = forward).
fn apply_strand_and_rc(seq: String, r: &Region, rc: bool) -> String {
    let strand_rc = r.strand == Some('-');
    // XOR: if both strand_rc and global rc are set, they cancel out.
    if strand_rc ^ rc {
        reverse_complement(&seq)
    } else {
        seq
    }
}

/// Apply masking and transforms to an extracted sequence.
fn apply_post_processing(seq: String, r: &Region, config: &OutputConfig) -> String {
    let mut s = seq;

    // Apply masking.
    if let Some(ref mask_idx) = config.mask_index {
        s = mask_idx.apply(&s, &r.chr, r.start, r.end, config.mask_mode);
    }

    // Apply transforms.
    if config.transform.any_active() {
        s = config.transform.apply(&s);
    }

    s
}

/// Resolve the correct FASTA path for a given region.
fn fasta_for_region<'a>(
    fasta_paths: &'a [PathBuf],
    contig_to_fasta: &HashMap<String, usize>,
    chr: &str,
) -> &'a Path {
    if fasta_paths.len() == 1 {
        return &fasta_paths[0];
    }
    let idx = contig_to_fasta.get(chr).copied().unwrap_or(0);
    &fasta_paths[idx]
}

/// Create a progress bar for extraction.
fn make_progress_bar(total: u64, show: bool) -> ProgressBar {
    if !show {
        return ProgressBar::hidden();
    }
    let pb = ProgressBar::new(total);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("[{elapsed_precise}] [{bar:40}] {pos}/{len} regions ({eta})")
            .expect("invalid progress bar template")
            .progress_chars("=>-"),
    );
    pb
}

#[allow(clippy::too_many_arguments)]
fn extract_and_format(
    fasta_paths: &[PathBuf],
    contig_to_fasta: &HashMap<String, usize>,
    fai_index: &HashMap<String, FaiRecord>,
    r: &Region,
    rc: bool,
    line_width: usize,
    index: usize,
    config: &OutputConfig,
) -> Result<(Region, String)> {
    let fasta_path = fasta_for_region(fasta_paths, contig_to_fasta, &r.chr);
    thread_local! {
        // Cache one file handle per FASTA path per thread.
        static TL_FILES: RefCell<HashMap<PathBuf, File>> = RefCell::new(HashMap::new());
    }
    TL_FILES.with(|cell| {
        let mut map = cell.borrow_mut();
        if !map.contains_key(fasta_path) {
            let f = File::open(fasta_path)
                .with_context(|| format!("Cannot open FASTA: {}", fasta_path.display()))?;
            map.insert(fasta_path.to_path_buf(), f);
        }
        let f = map.get_mut(fasta_path).unwrap();
        let seq = extract_region(f, fai_index, r)?;
        let seq = apply_strand_and_rc(seq, r, rc);
        let seq = apply_post_processing(seq, r, config);
        let entry = format_entry_with_config(r, &seq, line_width, index, config);
        Ok((r.clone(), entry))
    })
}

#[allow(clippy::too_many_arguments)]
pub fn write_sequences(
    fasta_paths: &[PathBuf],
    fai_index: &HashMap<String, FaiRecord>,
    contig_to_fasta: &HashMap<String, usize>,
    regions: &[Region],
    output_file: Option<&Path>,
    output_dir: Option<&Path>,
    is_sv: bool,
    rc: bool,
    line_width: usize,
    tab_out: bool,
    config: &OutputConfig,
    show_progress: bool,
) -> Result<()> {
    match output_dir {
        Some(dir) => {
            std::fs::create_dir_all(dir)?;
            if is_sv {
                write_sv_per_file_streaming(
                    fasta_paths,
                    contig_to_fasta,
                    fai_index,
                    regions,
                    dir,
                    rc,
                    line_width,
                    config,
                    show_progress,
                )?;
            } else {
                write_per_file_streaming(
                    fasta_paths,
                    contig_to_fasta,
                    fai_index,
                    regions,
                    dir,
                    rc,
                    line_width,
                    config,
                    show_progress,
                )?;
            }
        }
        None => {
            if tab_out {
                write_tab_output(
                    fasta_paths,
                    contig_to_fasta,
                    fai_index,
                    regions,
                    output_file,
                    rc,
                    config,
                    show_progress,
                )?;
            } else {
                write_streaming_ordered(
                    fasta_paths,
                    contig_to_fasta,
                    fai_index,
                    regions,
                    output_file,
                    rc,
                    line_width,
                    config,
                    show_progress,
                )?;
            }
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
#[allow(clippy::too_many_arguments)]
fn write_streaming_ordered(
    fasta_paths: &[PathBuf],
    contig_to_fasta: &HashMap<String, usize>,
    fai_index: &HashMap<String, FaiRecord>,
    regions: &[Region],
    output_file: Option<&Path>,
    rc: bool,
    line_width: usize,
    config: &OutputConfig,
    show_progress: bool,
) -> Result<()> {
    // When using a single FASTA file, use the optimised bulk-read path.
    // For multiple FASTA files, fall back to per-region extraction.
    if fasta_paths.len() == 1 {
        return write_streaming_ordered_single(
            &fasta_paths[0],
            fai_index,
            regions,
            output_file,
            rc,
            line_width,
            config,
            show_progress,
        );
    }

    let pb = make_progress_bar(regions.len() as u64, show_progress);
    let slots: Vec<Mutex<Option<String>>> = (0..regions.len()).map(|_| Mutex::new(None)).collect();

    regions
        .par_iter()
        .enumerate()
        .try_for_each(|(i, r)| -> Result<()> {
            let (_, entry) = extract_and_format(
                fasta_paths,
                contig_to_fasta,
                fai_index,
                r,
                rc,
                line_width,
                i,
                config,
            )?;
            let mut slot = slots[i].lock().unwrap();
            *slot = Some(entry);
            pb.inc(1);
            Ok(())
        })?;

    pb.finish_and_clear();

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

/// Optimised path for a single FASTA file using bulk reads.
#[allow(clippy::too_many_arguments)]
fn write_streaming_ordered_single(
    fasta_path: &Path,
    fai_index: &HashMap<String, FaiRecord>,
    regions: &[Region],
    output_file: Option<&Path>,
    rc: bool,
    line_width: usize,
    config: &OutputConfig,
    show_progress: bool,
) -> Result<()> {
    let groups = build_bulk_groups(regions, fai_index)?;
    let pb = make_progress_bar(regions.len() as u64, show_progress);
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
                let seq = apply_strand_and_rc(seq, r, rc);
                let seq = apply_post_processing(seq, r, config);
                let entry = format_entry_with_config(r, &seq, line_width, orig_idx, config);
                let mut slot = slots[orig_idx].lock().unwrap();
                *slot = Some(entry);
                pb.inc(1);
            }
            Ok(())
        })
    })?;
    pb.finish_and_clear();
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

/// Write TSV output: chr, start, end, name, sequence.
#[allow(clippy::too_many_arguments)]
fn write_tab_output(
    fasta_paths: &[PathBuf],
    contig_to_fasta: &HashMap<String, usize>,
    fai_index: &HashMap<String, FaiRecord>,
    regions: &[Region],
    output_file: Option<&Path>,
    rc: bool,
    config: &OutputConfig,
    show_progress: bool,
) -> Result<()> {
    if fasta_paths.len() == 1 {
        return write_tab_output_single(
            &fasta_paths[0],
            fai_index,
            regions,
            output_file,
            rc,
            config,
            show_progress,
        );
    }

    let pb = make_progress_bar(regions.len() as u64, show_progress);
    let slots: Vec<Mutex<Option<String>>> = (0..regions.len()).map(|_| Mutex::new(None)).collect();

    regions
        .par_iter()
        .enumerate()
        .try_for_each(|(i, r)| -> Result<()> {
            let fasta_path = fasta_for_region(fasta_paths, contig_to_fasta, &r.chr);
            thread_local! {
                static TL_FILES: RefCell<HashMap<PathBuf, File>> = RefCell::new(HashMap::new());
            }
            TL_FILES.with(|cell| {
                let mut map = cell.borrow_mut();
                if !map.contains_key(fasta_path) {
                    let f = File::open(fasta_path)
                        .with_context(|| format!("Cannot open FASTA: {}", fasta_path.display()))?;
                    map.insert(fasta_path.to_path_buf(), f);
                }
                let f = map.get_mut(fasta_path).unwrap();
                let seq = extract_region(f, fai_index, r)?;
                let seq = apply_strand_and_rc(seq, r, rc);
                let seq = apply_post_processing(seq, r, config);
                let entry = format_tab_entry(r, &seq);
                let mut slot = slots[i].lock().unwrap();
                *slot = Some(entry);
                pb.inc(1);
                Ok(())
            })
        })?;

    pb.finish_and_clear();

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

/// Optimised tab output path for a single FASTA file.
#[allow(clippy::too_many_arguments)]
fn write_tab_output_single(
    fasta_path: &Path,
    fai_index: &HashMap<String, FaiRecord>,
    regions: &[Region],
    output_file: Option<&Path>,
    rc: bool,
    config: &OutputConfig,
    show_progress: bool,
) -> Result<()> {
    let groups = build_bulk_groups(regions, fai_index)?;
    let pb = make_progress_bar(regions.len() as u64, show_progress);
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
                let seq = apply_strand_and_rc(seq, r, rc);
                let seq = apply_post_processing(seq, r, config);
                let entry = format_tab_entry(r, &seq);
                let mut slot = slots[orig_idx].lock().unwrap();
                *slot = Some(entry);
                pb.inc(1);
            }
            Ok(())
        })
    })?;
    pb.finish_and_clear();
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

#[allow(clippy::too_many_arguments)]
fn write_per_file_streaming(
    fasta_paths: &[PathBuf],
    contig_to_fasta: &HashMap<String, usize>,
    fai_index: &HashMap<String, FaiRecord>,
    regions: &[Region],
    dir: &Path,
    rc: bool,
    line_width: usize,
    config: &OutputConfig,
    show_progress: bool,
) -> Result<()> {
    let pb = make_progress_bar(regions.len() as u64, show_progress);
    regions
        .par_iter()
        .enumerate()
        .try_for_each(|(i, r)| -> Result<()> {
            let (region, entry) = extract_and_format(
                fasta_paths,
                contig_to_fasta,
                fai_index,
                r,
                rc,
                line_width,
                i,
                config,
            )?;
            let filename = match &region.name {
                Some(name) => format!("{name}_{}_{}.fa", region.start, region.end),
                None => format!("{}_{}_{}.fa", region.chr, region.start, region.end),
            };
            let mut f = File::create(dir.join(filename))?;
            f.write_all(entry.as_bytes())?;
            pb.inc(1);
            Ok(())
        })?;
    pb.finish_and_clear();
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn write_sv_per_file_streaming(
    fasta_paths: &[PathBuf],
    contig_to_fasta: &HashMap<String, usize>,
    fai_index: &HashMap<String, FaiRecord>,
    regions: &[Region],
    dir: &Path,
    rc: bool,
    line_width: usize,
    config: &OutputConfig,
    show_progress: bool,
) -> Result<()> {
    if !regions.len().is_multiple_of(2) {
        return Err(anyhow!(
            "SV extraction produced {} regions (expected even number for pairs)",
            regions.len()
        ));
    }
    let pb = make_progress_bar(regions.len() as u64, show_progress);
    regions
        .chunks(2)
        .enumerate()
        .collect::<Vec<_>>()
        .par_iter()
        .try_for_each(|&(chunk_idx, pair)| -> Result<()> {
            let idx1 = chunk_idx * 2;
            let idx2 = chunk_idx * 2 + 1;
            let (r1, entry1) = extract_and_format(
                fasta_paths,
                contig_to_fasta,
                fai_index,
                &pair[0],
                rc,
                line_width,
                idx1,
                config,
            )?;
            let (r2, entry2) = extract_and_format(
                fasta_paths,
                contig_to_fasta,
                fai_index,
                &pair[1],
                rc,
                line_width,
                idx2,
                config,
            )?;
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
            pb.inc(2);
            Ok(())
        })?;
    pb.finish_and_clear();
    Ok(())
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
    fn strand_suffix_plus() {
        let r = Region {
            name: None,
            chr: "chr1".into(),
            start: 1,
            end: 10,
            strand: Some('+'),
        };
        assert_eq!(strand_suffix(&r), "(+)");
    }

    #[test]
    fn strand_suffix_minus() {
        let r = Region {
            name: None,
            chr: "chr1".into(),
            start: 1,
            end: 10,
            strand: Some('-'),
        };
        assert_eq!(strand_suffix(&r), "(-)");
    }

    #[test]
    fn strand_suffix_none() {
        let r = Region {
            name: None,
            chr: "chr1".into(),
            start: 1,
            end: 10,
            strand: None,
        };
        assert_eq!(strand_suffix(&r), "");
    }

    #[test]
    fn apply_strand_and_rc_minus_strand_only() {
        // Minus strand alone should reverse complement.
        let r = Region {
            name: None,
            chr: "chr1".into(),
            start: 1,
            end: 4,
            strand: Some('-'),
        };
        assert_eq!(apply_strand_and_rc("ACGT".to_string(), &r, false), "ACGT");
        // ACGT RC = ACGT (palindrome). Use a non-palindrome.
        assert_eq!(apply_strand_and_rc("AAAC".to_string(), &r, false), "GTTT");
    }

    #[test]
    fn apply_strand_and_rc_cancel() {
        // Minus strand + global RC cancel out (XOR).
        let r = Region {
            name: None,
            chr: "chr1".into(),
            start: 1,
            end: 4,
            strand: Some('-'),
        };
        assert_eq!(apply_strand_and_rc("AAAC".to_string(), &r, true), "AAAC");
    }

    #[test]
    fn format_fasta_entry_with_strand() {
        let r = Region {
            name: Some("gene1".into()),
            chr: "chr1".into(),
            start: 1,
            end: 10,
            strand: Some('+'),
        };
        let entry = format_fasta_entry(&r, "AAACCCGGGT", 60);
        assert!(entry.starts_with(">gene1 chr1:1-10(+)\n"));
    }

    #[test]
    fn format_tab_entry_basic() {
        let r = Region {
            name: Some("gene1".into()),
            chr: "chr1".into(),
            start: 1,
            end: 10,
            strand: None,
        };
        let entry = format_tab_entry(&r, "AAACCCGGGT");
        assert_eq!(entry, "chr1\t1\t10\tgene1\tAAACCCGGGT\n");
    }

    #[test]
    fn format_tab_entry_no_name() {
        let r = Region {
            name: None,
            chr: "chr1".into(),
            start: 1,
            end: 10,
            strand: None,
        };
        let entry = format_tab_entry(&r, "AAACCCGGGT");
        assert_eq!(entry, "chr1\t1\t10\t.\tAAACCCGGGT\n");
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
            strand: None,
        })
        .collect();
        let out_file = tempfile::NamedTempFile::new().unwrap();
        let config = OutputConfig::default();
        write_streaming_ordered_single(
            tmp.path(),
            &fai,
            &regions,
            Some(out_file.path()),
            false,
            60,
            &config,
            false,
        )
        .unwrap();
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
