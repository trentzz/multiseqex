use anyhow::{Context, Result, anyhow};
use clap::Parser;
use rayon::prelude::*;
use std::cmp::{max, min};
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};

// ─── Data types ──────────────────────────────────────────────────────────────

/// A genomic interval (1-based, inclusive on both ends).
#[derive(Debug, Clone)]
struct Region {
    name: Option<String>,
    chr: String,
    start: u64,
    end: u64,
}

/// One record from a `.fai` index file.
#[derive(Debug, Clone)]
struct FaiRecord {
    /// Total number of bases in this contig.
    length: u64,
    /// Byte offset of the first base in the FASTA file.
    offset: u64,
    /// Number of sequence bases per line.
    line_bases: u64,
    /// Number of bytes per line (bases + newline characters).
    line_bytes: u64,
}

// ─── CLI definition ──────────────────────────────────────────────────────────

#[derive(Parser, Debug)]
#[command(
    name = "multiseqex",
    author,
    version,
    about = "Multi-sequence extractor for FASTA using FAI"
)]
struct Cli {
    /// Reference FASTA file (bgzipped ok if a matching .fai exists).
    fasta: PathBuf,

    /// Comma-separated regions: chr:start-end, chr2:start-end, ...
    #[arg(long)]
    regions: Option<String>,

    /// File with one region per line (chr:start-end).
    #[arg(long)]
    list: Option<PathBuf>,

    /// CSV/TSV table with named columns.
    ///
    /// Required: CHROM.
    /// Range mode: CHROM, START, END.
    /// Position mode: CHROM, POS (requires --flank).
    /// Optional: NAME.
    /// Extra columns are ignored. Delimiter auto-detected (.tsv → tab, else comma).
    #[arg(long)]
    table: Option<PathBuf>,

    /// CSV/TSV structural-variant table with named columns.
    ///
    /// Required: CHROM_LEFT, CHROM_RIGHT.
    /// Range mode: START_LEFT, END_LEFT, START_RIGHT, END_RIGHT.
    /// Position mode: POS_LEFT, POS_RIGHT (requires --flank).
    /// Optional: NAME.
    /// Extra columns are ignored.
    #[arg(long)]
    sv_table: Option<PathBuf>,

    /// Flank size for position-mode tables (required with POS columns).
    #[arg(long)]
    flank: Option<u64>,

    /// Output FASTA file (single combined file; default: stdout).
    #[arg(short, long)]
    output: Option<PathBuf>,

    /// Output directory — one FASTA file per region (or per SV pair).
    #[arg(long)]
    output_dir: Option<PathBuf>,

    /// Number of worker threads (default: all available CPUs).
    #[arg(long)]
    threads: Option<usize>,

    /// Error if .fai is missing instead of building one automatically.
    #[arg(long)]
    no_build_fai: bool,
}

// ─── Entry point ─────────────────────────────────────────────────────────────

fn main() -> Result<()> {
    let cli = Cli::parse();

    if let Some(t) = cli.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(t)
            .build_global()
            .ok();
    }

    if cli.output.is_some() && cli.output_dir.is_some() {
        return Err(anyhow!(
            "Cannot use both --output and --output-dir simultaneously"
        ));
    }

    // Ensure .fai exists (or build it).
    let fai_path = fai_path_for(&cli.fasta);
    if !fai_path.exists() {
        if cli.no_build_fai {
            return Err(anyhow!(
                "Missing index: {} (use samtools faidx or remove --no-build-fai)",
                fai_path.display()
            ));
        }
        eprintln!("Index not found. Building FAI: {}", fai_path.display());
        build_fai(&cli.fasta, &fai_path)?;
    }
    let fai_index = read_fai(&fai_path)?;

    // Collect regions from all input sources.
    let mut regions = Vec::<Region>::new();
    if let Some(s) = cli.regions.as_deref() {
        regions.extend(parse_regions_inline(s, cli.flank)?);
    }
    if let Some(p) = cli.list.as_ref() {
        regions.extend(parse_regions_list(p, cli.flank)?);
    }
    if let Some(p) = cli.table.as_ref() {
        regions.extend(parse_regions_table(p, cli.flank)?);
    }
    if let Some(p) = cli.sv_table.as_ref() {
        regions.extend(parse_regions_sv_table(p, cli.flank)?);
    }

    if regions.is_empty() {
        return Err(anyhow!(
            "No regions provided. Use --regions, --list, --table, or --sv-table."
        ));
    }

    // Validate and clamp regions to contig lengths.
    validate_and_clamp_regions(&mut regions, &fai_index)?;

    let is_sv = cli.sv_table.is_some();
    write_sequences(
        &cli.fasta,
        &fai_index,
        &regions,
        cli.output.as_deref(),
        cli.output_dir.as_deref(),
        is_sv,
    )?;

    Ok(())
}

// ─── FAI helpers ─────────────────────────────────────────────────────────────

/// Compute the `.fai` path for a given FASTA file.
fn fai_path_for(fasta: &Path) -> PathBuf {
    let mut s = fasta.as_os_str().to_owned();
    s.push(".fai");
    PathBuf::from(s)
}

/// Build a minimal `.fai` index from a FASTA file.
fn build_fai(fasta: &Path, fai_out: &Path) -> Result<()> {
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
            // EOF — flush last contig.
            if let Some(name) = current_name.take() {
                writeln!(out, "{name}\t{seq_len}\t{seq_offset}\t{line_bases}\t{line_bytes}")?;
            }
            break;
        }

        let raw = line.as_bytes();
        let linelen = raw.len() as u64;

        if raw.starts_with(b">") {
            // New contig header — flush the previous one.
            if let Some(name) = current_name.replace(parse_fasta_header(&line)) {
                writeln!(out, "{name}\t{seq_len}\t{seq_offset}\t{line_bases}\t{line_bytes}")?;
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
            }
        }
        pos += linelen;
    }

    out.flush()?;
    Ok(())
}

/// Extract the first whitespace-delimited token after `>`.
fn parse_fasta_header(s: &str) -> String {
    s.trim_start_matches('>')
        .split_whitespace()
        .next()
        .unwrap_or("")
        .to_string()
}

/// Count IUPAC nucleotide characters in a raw line.
fn count_bases(raw: &[u8]) -> u64 {
    raw.iter()
        .filter(|b| b.is_ascii_alphabetic())
        .count() as u64
}

/// Read a `.fai` file into a contig-name → `FaiRecord` map.
fn read_fai(fai_path: &Path) -> Result<HashMap<String, FaiRecord>> {
    let f = File::open(fai_path)
        .with_context(|| format!("Cannot open FAI: {}", fai_path.display()))?;
    let reader = BufReader::new(f);
    let mut map = HashMap::new();

    for (i, line) in reader.lines().enumerate() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 5 {
            return Err(anyhow!("Malformed FAI line {}: {}", i + 1, line));
        }
        map.insert(
            parts[0].to_string(),
            FaiRecord {
                length: parts[1].parse().with_context(|| format!("Bad length, FAI line {}", i + 1))?,
                offset: parts[2].parse().with_context(|| format!("Bad offset, FAI line {}", i + 1))?,
                line_bases: parts[3].parse().with_context(|| format!("Bad line_bases, FAI line {}", i + 1))?,
                line_bytes: parts[4].parse().with_context(|| format!("Bad line_bytes, FAI line {}", i + 1))?,
            },
        );
    }
    Ok(map)
}

// ─── Sequence extraction ─────────────────────────────────────────────────────

/// Extract a region from a FASTA file using the FAI index.
fn extract_region(fasta: &Path, fai: &HashMap<String, FaiRecord>, r: &Region) -> Result<String> {
    let rec = fai
        .get(&r.chr)
        .ok_or_else(|| anyhow!("Contig '{}' not in index", r.chr))?;

    if r.start < 1 || r.end < 1 || r.start > rec.length || r.end > rec.length {
        return Err(anyhow!(
            "Region out of bounds {}:{}-{} (contig length={})",
            r.chr,
            r.start,
            r.end,
            rec.length
        ));
    }

    let mut f =
        File::open(fasta).with_context(|| format!("Cannot open FASTA: {}", fasta.display()))?;
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

// ─── Region parsing: inline & list ───────────────────────────────────────────

/// Parse comma-separated inline region strings.
fn parse_regions_inline(s: &str, flank: Option<u64>) -> Result<Vec<Region>> {
    s.split(',')
        .filter(|t| !t.trim().is_empty())
        .map(|t| parse_region_str(t.trim(), flank))
        .collect()
}

/// Parse a file with one region string per line.
fn parse_regions_list(path: &Path, flank: Option<u64>) -> Result<Vec<Region>> {
    let f = File::open(path)
        .with_context(|| format!("Cannot open list file: {}", path.display()))?;
    BufReader::new(f)
        .lines()
        .filter_map(|l| {
            let l = l.ok()?;
            let trimmed = l.trim().to_string();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                None
            } else {
                Some(trimmed)
            }
        })
        .map(|l| parse_region_str(&l, flank))
        .collect()
}

/// Parse a single region string: `chr:start-end` or `chr:pos+flank`.
fn parse_region_str(s: &str, flank: Option<u64>) -> Result<Region> {
    let (chr, rest) = s
        .split_once(':')
        .ok_or_else(|| anyhow!("Bad region (missing ':'): {}", s))?;

    // chr:start-end
    if let Some((start_s, end_s)) = rest.split_once('-') {
        let start: u64 = start_s.replace(',', "").parse()
            .with_context(|| format!("Bad start in region: {s}"))?;
        let end: u64 = end_s.replace(',', "").parse()
            .with_context(|| format!("Bad end in region: {s}"))?;
        return Ok(Region {
            name: None,
            chr: chr.to_string(),
            start: min(start, end),
            end: max(start, end),
        });
    }

    // chr:pos+flank
    let (pos_s, flank_s) = rest
        .split_once('+')
        .ok_or_else(|| anyhow!("Bad region format (expected start-end or pos+flank): {s}"))?;

    let pos: u64 = pos_s.replace(',', "").parse()
        .with_context(|| format!("Bad position in region: {s}"))?;
    let inline_flank: u64 = flank_s.replace(',', "").parse()
        .with_context(|| format!("Bad flank in region: {s}"))?;
    let effective_flank = inline_flank.max(flank.unwrap_or(0));

    Ok(Region {
        name: None,
        chr: chr.to_string(),
        start: pos.saturating_sub(effective_flank).max(1),
        end: pos + effective_flank,
    })
}

// ─── Region parsing: CSV/TSV tables ──────────────────────────────────────────

/// Build a case-insensitive header-name → column-index map.
fn build_header_map(headers: &csv::StringRecord) -> HashMap<String, usize> {
    headers
        .iter()
        .enumerate()
        .map(|(i, h)| (h.trim().to_uppercase(), i))
        .collect()
}

/// Auto-detect CSV/TSV delimiter from file extension.
fn detect_delimiter(path: &Path) -> u8 {
    let is_tsv = path
        .extension()
        .and_then(|e| e.to_str())
        .is_some_and(|e| e.eq_ignore_ascii_case("tsv"));
    if is_tsv { b'\t' } else { b',' }
}

/// Parse a numeric field from a CSV record, with a contextual error message.
fn parse_u64_field(
    rec: &csv::StringRecord,
    col_idx: usize,
    field_name: &str,
    row: usize,
) -> Result<u64> {
    rec.get(col_idx)
        .ok_or_else(|| anyhow!("Missing {field_name} at row {row}"))?
        .trim()
        .replace(',', "")
        .parse()
        .with_context(|| format!("Bad {field_name} at row {row}"))
}

/// Look up a required column by name, returning a helpful error if absent.
fn require_column(
    hmap: &HashMap<String, usize>,
    name: &str,
    headers: &csv::StringRecord,
) -> Result<usize> {
    hmap.get(name)
        .copied()
        .ok_or_else(|| anyhow!("Table missing required {name} column (found: {headers:?})"))
}

/// Read a string field from a CSV record.
fn read_string_field(rec: &csv::StringRecord, col_idx: usize) -> Result<String> {
    Ok(rec
        .get(col_idx)
        .unwrap_or("")
        .trim()
        .to_string())
}

/// Read an optional NAME field (returns `None` if the column is absent or empty).
fn read_optional_name(rec: &csv::StringRecord, name_idx: Option<usize>) -> Option<String> {
    name_idx.and_then(|idx| {
        rec.get(idx)
            .map(|v| v.trim().to_string())
            .filter(|v| !v.is_empty())
    })
}

enum TableMode {
    Range { start_idx: usize, end_idx: usize },
    Position { pos_idx: usize },
}

/// Parse a CSV/TSV table with named columns: CHROM, START, END, POS, NAME.
fn parse_regions_table(path: &Path, flank: Option<u64>) -> Result<Vec<Region>> {
    let delim = detect_delimiter(path);
    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(true)
        .delimiter(delim)
        .from_path(path)
        .with_context(|| format!("Cannot open table: {}", path.display()))?;

    let headers = rdr.headers()?.clone();
    let hmap = build_header_map(&headers);

    let chrom_idx = require_column(&hmap, "CHROM", &headers)?;
    let name_idx = hmap.get("NAME").copied();

    let has_start = hmap.get("START").copied();
    let has_end = hmap.get("END").copied();
    let has_pos = hmap.get("POS").copied();

    let mode = match (has_start, has_end, has_pos) {
        (Some(s), Some(e), None) => TableMode::Range {
            start_idx: s,
            end_idx: e,
        },
        (None, None, Some(p)) => TableMode::Position { pos_idx: p },
        (Some(_), Some(_), Some(_)) => {
            return Err(anyhow!(
                "Table has both START/END and POS columns — ambiguous. \
                 Use either START+END (range) or POS (position)."
            ));
        }
        (Some(_), None, _) | (None, Some(_), _) => {
            return Err(anyhow!(
                "Table has only one of START/END — both are required for range mode (found: {headers:?})"
            ));
        }
        _ => {
            return Err(anyhow!(
                "Table must have START+END or POS columns (found: {headers:?})"
            ));
        }
    };

    if matches!(mode, TableMode::Position { .. }) && flank.is_none() {
        return Err(anyhow!(
            "--flank is required when table uses POS column (position mode)"
        ));
    }
    let flank = flank.unwrap_or(0);

    let mut out = Vec::new();
    for (i, rec) in rdr.records().enumerate() {
        let rec = rec?;
        let row = i + 2;
        let chr = read_string_field(&rec, chrom_idx)?;
        let name = read_optional_name(&rec, name_idx);

        let region = match &mode {
            TableMode::Range { start_idx, end_idx } => {
                let s = parse_u64_field(&rec, *start_idx, "START", row)?;
                let e = parse_u64_field(&rec, *end_idx, "END", row)?;
                Region { name, chr, start: min(s, e), end: max(s, e) }
            }
            TableMode::Position { pos_idx } => {
                let p = parse_u64_field(&rec, *pos_idx, "POS", row)?;
                Region { name, chr, start: p.saturating_sub(flank).max(1), end: p + flank }
            }
        };
        out.push(region);
    }
    Ok(out)
}

enum SvMode {
    Range {
        start_left_idx: usize,
        end_left_idx: usize,
        start_right_idx: usize,
        end_right_idx: usize,
    },
    Position {
        pos_left_idx: usize,
        pos_right_idx: usize,
    },
}

/// Parse a CSV/TSV SV table with named columns.
fn parse_regions_sv_table(path: &Path, flank: Option<u64>) -> Result<Vec<Region>> {
    let delim = detect_delimiter(path);
    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(true)
        .delimiter(delim)
        .from_path(path)
        .with_context(|| format!("Cannot open SV table: {}", path.display()))?;

    let headers = rdr.headers()?.clone();
    let hmap = build_header_map(&headers);

    let chrom_left_idx = require_column(&hmap, "CHROM_LEFT", &headers)?;
    let chrom_right_idx = require_column(&hmap, "CHROM_RIGHT", &headers)?;
    let name_idx = hmap.get("NAME").copied();

    let has_sl = hmap.get("START_LEFT").copied();
    let has_el = hmap.get("END_LEFT").copied();
    let has_sr = hmap.get("START_RIGHT").copied();
    let has_er = hmap.get("END_RIGHT").copied();
    let has_pl = hmap.get("POS_LEFT").copied();
    let has_pr = hmap.get("POS_RIGHT").copied();

    let mode = match (has_sl, has_el, has_sr, has_er, has_pl, has_pr) {
        (Some(sl), Some(el), Some(sr), Some(er), _, _) => SvMode::Range {
            start_left_idx: sl,
            end_left_idx: el,
            start_right_idx: sr,
            end_right_idx: er,
        },
        (None, None, None, None, Some(pl), Some(pr)) => SvMode::Position {
            pos_left_idx: pl,
            pos_right_idx: pr,
        },
        _ => {
            return Err(anyhow!(
                "SV table must have START_LEFT+END_LEFT+START_RIGHT+END_RIGHT (range) \
                 or POS_LEFT+POS_RIGHT (position). Found: {headers:?}"
            ));
        }
    };

    if matches!(mode, SvMode::Position { .. }) && flank.is_none() {
        return Err(anyhow!(
            "--flank is required when SV table uses POS_LEFT/POS_RIGHT (position mode)"
        ));
    }
    let flank = flank.unwrap_or(0);

    let mut out = Vec::new();
    for (i, rec) in rdr.records().enumerate() {
        let rec = rec?;
        let row = i + 2;
        let chr_left = read_string_field(&rec, chrom_left_idx)?;
        let chr_right = read_string_field(&rec, chrom_right_idx)?;
        let name = read_optional_name(&rec, name_idx);

        match &mode {
            SvMode::Range { start_left_idx, end_left_idx, start_right_idx, end_right_idx } => {
                let sl = parse_u64_field(&rec, *start_left_idx, "START_LEFT", row)?;
                let el = parse_u64_field(&rec, *end_left_idx, "END_LEFT", row)?;
                let sr = parse_u64_field(&rec, *start_right_idx, "START_RIGHT", row)?;
                let er = parse_u64_field(&rec, *end_right_idx, "END_RIGHT", row)?;
                out.push(Region { name: name.clone(), chr: chr_left, start: min(sl, el), end: max(sl, el) });
                out.push(Region { name, chr: chr_right, start: min(sr, er), end: max(sr, er) });
            }
            SvMode::Position { pos_left_idx, pos_right_idx } => {
                let pl = parse_u64_field(&rec, *pos_left_idx, "POS_LEFT", row)?;
                let pr = parse_u64_field(&rec, *pos_right_idx, "POS_RIGHT", row)?;
                out.push(Region { name: name.clone(), chr: chr_left, start: pl.saturating_sub(flank).max(1), end: pl + flank });
                out.push(Region { name, chr: chr_right, start: pr.saturating_sub(flank).max(1), end: pr + flank });
            }
        }
    }
    Ok(out)
}

// ─── Validation ──────────────────────────────────────────────────────────────

/// Validate that all regions reference known contigs and clamp to contig bounds.
fn validate_and_clamp_regions(
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

// ─── Output ──────────────────────────────────────────────────────────────────

/// Wrap a sequence string to a fixed line width for FASTA output.
fn wrap_fasta(seq: &str, width: usize) -> String {
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

/// Extract and write all sequences to the requested destination.
fn write_sequences(
    fasta_path: &Path,
    fai_index: &HashMap<String, FaiRecord>,
    regions: &[Region],
    output_file: Option<&Path>,
    output_dir: Option<&Path>,
    is_sv: bool,
) -> Result<()> {
    // Parallel extraction.
    let sequences: Vec<(Region, String)> = regions
        .par_iter()
        .map(|r| {
            let seq = extract_region(fasta_path, fai_index, r)?;
            let header = match &r.name {
                Some(name) => format!(">{name} {}:{}-{}", r.chr, r.start, r.end),
                None => format!(">{}:{}-{}", r.chr, r.start, r.end),
            };
            let fasta_entry = format!(
                "{header}\n{}\n",
                wrap_fasta(&seq, 60)
            );
            Ok((r.clone(), fasta_entry))
        })
        .collect::<Result<_>>()?;

    match output_dir {
        Some(dir) => {
            std::fs::create_dir_all(dir)?;
            if is_sv {
                write_sv_per_file(dir, &sequences)?;
            } else {
                write_per_file(dir, &sequences)?;
            }
        }
        None => {
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

/// Write one FASTA file per region.
fn write_per_file(dir: &Path, sequences: &[(Region, String)]) -> Result<()> {
    sequences
        .par_iter()
        .try_for_each(|(region, entry)| -> Result<()> {
            let filename = match &region.name {
                Some(name) => format!("{name}_{}_{}.fa", region.start, region.end),
                None => format!("{}_{}_{}.fa", region.chr, region.start, region.end),
            };
            let mut f = File::create(dir.join(filename))?;
            f.write_all(entry.as_bytes())?;
            Ok(())
        })
}

/// Write one FASTA file per SV pair (two regions per file).
fn write_sv_per_file(dir: &Path, sequences: &[(Region, String)]) -> Result<()> {
    if !sequences.len().is_multiple_of(2) {
        return Err(anyhow!(
            "SV extraction produced {} regions (expected even number for pairs)",
            sequences.len()
        ));
    }

    sequences
        .chunks(2)
        .collect::<Vec<_>>()
        .par_iter()
        .try_for_each(|pair| -> Result<()> {
            let (r1, _) = &pair[0];
            let (r2, _) = &pair[1];

            let filename = match (&r1.name, &r2.name) {
                (Some(n1), Some(n2)) if n1 == n2 => {
                    format!("{n1}_{}_{}_{}_{}_{}_{}.fa", r1.chr, r1.start, r1.end, r2.chr, r2.start, r2.end)
                }
                (Some(n1), Some(n2)) => {
                    return Err(anyhow!("Mismatched names in SV pair: '{n1}' vs '{n2}'"));
                }
                _ => {
                    format!("{}_{}_{}_{}_{}_{}.fa", r1.chr, r1.start, r1.end, r2.chr, r2.start, r2.end)
                }
            };

            let mut f = File::create(dir.join(filename))?;
            for (_, entry) in pair.iter() {
                f.write_all(entry.as_bytes())?;
            }
            Ok(())
        })
}
