use anyhow::{Context, Result, anyhow};
use clap::{ArgAction, Parser};
use rayon::prelude::*;
use std::cmp::{max, min};
use std::collections::HashMap;
use std::fs::File;
use std::io::Read;
use std::io::Seek;
use std::io::{self, BufRead, BufReader, BufWriter, SeekFrom, Write};
use std::path::{Path, PathBuf};

/// A genomic interval (1-based inclusive).
#[derive(Debug, Clone)]
struct Region {
    name: Option<String>,
    chr: String,
    start: u64, // 1-based inclusive
    end: u64,   // 1-based inclusive
}

#[derive(Debug, Clone)]
struct FaiRecord {
    length: u64,
    offset: u64,     // byte offset of first base of this contig in the FASTA
    line_bases: u64, // number of bases per FASTA line
    line_bytes: u64, // number of bytes per FASTA line (includes newline(s))
}

/// CLI definition
#[derive(Parser, Debug)]
#[command(
    name = "multiseqex",
    author,
    version,
    about = "Multi-sequence extractor for FASTA using FAI"
)]
struct Cli {
    /// Reference FASTA file (bgzipped FASTA is fine as long as a matching .fai exists)
    fasta: PathBuf,

    /// Comma-separated regions: chr:start-end,chr2:start-end,...
    #[arg(long)]
    regions: Option<String>,

    /// File with one region per line in the form chr:start-end
    #[arg(long)]
    list: Option<PathBuf>,

    /// CSV/TSV table with named columns: CHROM + START/END (range mode) or CHROM + POS (position mode, requires --flank). Optional: NAME. Extra columns are ignored.
    /// Delimiter auto-detected from extension: .tsv => tab, else comma
    #[arg(long)]
    table: Option<PathBuf>,

    /// Flank size for position-mode tables (required when using POS or POS_LEFT/POS_RIGHT columns)
    #[arg(long)]
    flank: Option<u64>,

    /// Output FASTA file (single combined; default: stdout)
    #[arg(short, long)]
    output: Option<PathBuf>,

    /// Output directory for per-sequence FASTA files
    #[arg(long)]
    output_dir: Option<PathBuf>,

    /// Number of worker threads (default: Rayon default, usually #cpus)
    #[arg(long)]
    threads: Option<usize>,

    /// CSV/TSV SV table with named columns: CHROM_LEFT + CHROM_RIGHT + START_LEFT/END_LEFT/START_RIGHT/END_RIGHT (range mode) or POS_LEFT/POS_RIGHT (position mode, requires --flank). Optional: NAME. Extra columns are ignored.
    #[arg(long)]
    sv_table: Option<PathBuf>,

    /// Do not build a .fai if missing (error instead)
    #[arg(long, action=ArgAction::SetTrue)]
    no_build_fai: bool,
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    if let Some(t) = cli.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(t)
            .build_global()
            .ok();
    }

    // Ensure .fai exists (or build it)
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

    // Collect regions from any combination of inputs
    let mut regions: Vec<Region> = Vec::new();
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
            "No regions provided. Use --regions, --list, or --table."
        ));
    }

    // Validate & clamp to contig lengths
    let mut bad: Vec<String> = Vec::new();
    for r in regions.iter_mut() {
        if let Some(rec) = fai_index.get(&r.chr) {
            if r.start == 0 {
                r.start = 1;
            }
            if r.end == 0 {
                r.end = 1;
            }
            r.start = min(r.start, rec.length);
            r.end = min(r.end, rec.length);
            if r.start > r.end {
                std::mem::swap(&mut r.start, &mut r.end);
            }
        } else {
            bad.push(r.chr.clone());
        }
    }
    if !bad.is_empty() {
        return Err(anyhow!(
            "Contigs not found in FASTA/FAI: {}",
            unique_join(&bad, ", ")
        ));
    }

    // Enforce mutually exclusive output options
    if cli.output.is_some() && cli.output_dir.is_some() {
        return Err(anyhow!(
            "Cannot use both --output and --output-dir simultaneously"
        ));
    }

    // Extract and write sequences
    write_sequences(
        &cli.fasta,
        &fai_index,
        &regions,
        cli.output.as_deref(),
        cli.output_dir.as_deref(),
        cli.sv_table.as_deref(),
    )?;

    Ok(())
}

/// Compute the .fai path next to a FASTA
fn fai_path_for(fasta: &Path) -> PathBuf {
    let s = fasta.to_string_lossy();
    PathBuf::from(format!("{}.fai", s))
}

/// Minimal .fai builder
fn build_fai(fasta: &Path, fai_out: &Path) -> Result<()> {
    let f = File::open(fasta)
        .with_context(|| format!("Cannot open FASTA for indexing: {}", fasta.display()))?;
    let mut r = BufReader::new(f);
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
    let mut first_seq_line_seen = false;

    loop {
        line.clear();
        let n = r.read_line(&mut line)?;
        if n == 0 {
            if let Some(name) = current_name.take() {
                writeln!(
                    out,
                    "{}\t{}\t{}\t{}\t{}",
                    name, seq_len, seq_offset, line_bases, line_bytes
                )?;
            }
            break;
        }

        let raw = line.as_bytes();
        let linelen = raw.len() as u64;
        if raw.starts_with(b">") {
            if let Some(name) = current_name.replace(parse_fasta_header(&line)) {
                writeln!(
                    out,
                    "{}\t{}\t{}\t{}\t{}",
                    name, seq_len, seq_offset, line_bases, line_bytes
                )?;
            }
            seq_offset = pos + linelen;
            seq_len = 0;
            line_bases = 0;
            line_bytes = 0;
            first_seq_line_seen = false;
        } else {
            let (bases_count, bytes_count) = count_bases_and_bytes(raw);
            seq_len += bases_count;
            if !first_seq_line_seen {
                line_bases = bases_count;
                line_bytes = bytes_count;
                first_seq_line_seen = true;
            }
        }
        pos += linelen;
    }

    out.flush()?;
    Ok(())
}

fn parse_fasta_header(s: &str) -> String {
    s.trim_start_matches('>')
        .trim()
        .split_whitespace()
        .next()
        .unwrap_or("")
        .to_string()
}

fn count_bases_and_bytes(raw: &[u8]) -> (u64, u64) {
    let mut bases = 0u64;
    let bytes = raw.len() as u64;
    for &b in raw {
        if matches!(
            b,
            b'A' | b'C'
                | b'G'
                | b'T'
                | b'U'
                | b'R'
                | b'Y'
                | b'S'
                | b'W'
                | b'K'
                | b'M'
                | b'B'
                | b'D'
                | b'H'
                | b'V'
                | b'N'
                | b'a'
                | b'c'
                | b'g'
                | b't'
                | b'u'
                | b'r'
                | b'y'
                | b's'
                | b'w'
                | b'k'
                | b'm'
                | b'b'
                | b'd'
                | b'h'
                | b'v'
                | b'n'
        ) {
            bases += 1;
        }
    }
    (bases, bytes)
}

fn read_fai(fai_path: &Path) -> Result<HashMap<String, FaiRecord>> {
    let f =
        File::open(fai_path).with_context(|| format!("Cannot open FAI: {}", fai_path.display()))?;
    let r = BufReader::new(f);
    let mut map = HashMap::new();

    for (i, line) in r.lines().enumerate() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 5 {
            return Err(anyhow!("Malformed FAI line {}: {}", i + 1, line));
        }
        let name = parts[0].to_string();
        let length: u64 = parts[1].parse()?;
        let offset: u64 = parts[2].parse()?;
        let line_bases: u64 = parts[3].parse()?;
        let line_bytes: u64 = parts[4].parse()?;
        map.insert(
            name,
            FaiRecord {
                length,
                offset,
                line_bases,
                line_bytes,
            },
        );
    }
    Ok(map)
}

/// Extract a region using FAI math
fn extract_region(fasta: &Path, fai: &HashMap<String, FaiRecord>, r: &Region) -> Result<String> {
    let rec = fai
        .get(&r.chr)
        .ok_or_else(|| anyhow!("Contig {} not in index", r.chr))?;
    let start = r.start;
    let end = r.end;
    if start < 1 || end < 1 || start > rec.length || end > rec.length {
        return Err(anyhow!(
            "Region out of bounds {}:{}-{} (len={})",
            r.chr,
            start,
            end,
            rec.length
        ));
    }

    let mut f =
        File::open(fasta).with_context(|| format!("Cannot open FASTA: {}", fasta.display()))?;
    let lb = rec.line_bases;
    let lby = rec.line_bytes;

    let mut seq = Vec::<u8>::with_capacity((end - start + 1) as usize);
    let mut p = start;
    while p <= end {
        let line_idx = (p - 1) / lb;
        let in_line_offset = (p - 1) % lb;
        let to_line_end = lb - in_line_offset;
        let remaining = end - p + 1;
        let run = min(to_line_end, remaining);
        let byte_pos = rec.offset + line_idx * lby + in_line_offset;
        f.seek(SeekFrom::Start(byte_pos))?;
        let mut buf = vec![0u8; run as usize];
        f.read_exact(&mut buf)?;
        seq.extend_from_slice(&buf);
        p += run;
    }

    for b in seq.iter_mut() {
        *b = b.to_ascii_uppercase();
    }
    Ok(String::from_utf8(seq)?)
}

/// Parse comma-separated inline regions
fn parse_regions_inline(s: &str, flank: Option<u64>) -> Result<Vec<Region>> {
    let mut out = Vec::new();
    for (_i, tok) in s.split(',').enumerate() {
        if tok.trim().is_empty() {
            continue;
        }
        out.push(parse_region_str(tok.trim(), flank)?);
    }
    Ok(out)
}

/// Parse line-based list file
fn parse_regions_list(path: &Path, flank: Option<u64>) -> Result<Vec<Region>> {
    let f =
        File::open(path).with_context(|| format!("Cannot open list file: {}", path.display()))?;
    let r = BufReader::new(f);
    let mut out = Vec::new();
    for (_i, line) in r.lines().enumerate() {
        let l = line?;
        if l.trim().is_empty() || l.starts_with('#') {
            continue;
        }
        out.push(parse_region_str(l.trim(), flank)?);
    }
    Ok(out)
}

/// Parse a single region string
fn parse_region_str(s: &str, flank: Option<u64>) -> Result<Region> {
    // Case: chr:start-end or chr:pos±flank
    let (chr, rest) = s
        .split_once(':')
        .ok_or_else(|| anyhow!("Bad region (missing ':'): {}", s))?;
    if let Some((start, end)) = rest.split_once('-') {
        let start: u64 = start.replace(',', "").parse()?;
        let end: u64 = end.replace(',', "").parse()?;
        let (start, end) = if start <= end {
            (start, end)
        } else {
            (end, start)
        };
        return Ok(Region {
            name: None,
            chr: chr.to_string(),
            start,
            end,
        });
    }

    // pos±flank
    let (pos, flank_val) = if let Some((p, f)) = rest.split_once('+') {
        (p, f)
    } else if let Some((p, f)) = rest.split_once('-') {
        (p, f)
    } else {
        return Err(anyhow!(
            "Bad region format (expected start-end or pos±flank): {}",
            s
        ));
    };

    let flank = flank.unwrap_or(0);

    let pos: u64 = pos.replace(',', "").parse()?;
    let flank_val: u64 = flank_val.replace(',', "").parse()?;
    let start = pos.saturating_sub(flank_val.max(flank));
    let end = pos + flank_val.max(flank);

    Ok(Region {
        name: None,
        chr: chr.to_string(),
        start: start.max(1),
        end,
    })
}

/// Build a case-insensitive header name -> column index map
fn build_header_map(headers: &csv::StringRecord) -> HashMap<String, usize> {
    let mut map = HashMap::new();
    for (i, h) in headers.iter().enumerate() {
        map.insert(h.trim().to_uppercase(), i);
    }
    map
}

/// Auto-detect delimiter from file extension (.tsv -> tab, else comma)
fn detect_delimiter(path: &Path) -> u8 {
    if path
        .extension()
        .and_then(|e| e.to_str())
        .map(|e| e.eq_ignore_ascii_case("tsv"))
        .unwrap_or(false)
    {
        b'\t'
    } else {
        b','
    }
}

/// Parse a numeric field from a record by column index, with context on error
fn parse_u64_field(
    rec: &csv::StringRecord,
    col_idx: usize,
    field_name: &str,
    row: usize,
) -> Result<u64> {
    rec.get(col_idx)
        .ok_or_else(|| anyhow!("Missing {} at row {}", field_name, row))?
        .trim()
        .replace(',', "")
        .parse()
        .with_context(|| format!("Bad {} at row {}", field_name, row))
}

/// Parse CSV/TSV table using column names: CHROM, START, END, POS, NAME
fn parse_regions_table(path: &Path, flank: Option<u64>) -> Result<Vec<Region>> {
    let delim = detect_delimiter(path);
    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(true)
        .delimiter(delim)
        .from_path(path)
        .with_context(|| format!("Cannot open table: {}", path.display()))?;

    let headers = rdr.headers()?.clone();
    let hmap = build_header_map(&headers);

    // Required column
    let chrom_idx = hmap
        .get("CHROM")
        .copied()
        .ok_or_else(|| anyhow!("Table missing required CHROM column (found: {:?})", headers))?;

    // Detect mode from available columns
    let has_start = hmap.get("START").copied();
    let has_end = hmap.get("END").copied();
    let has_pos = hmap.get("POS").copied();
    let name_idx = hmap.get("NAME").copied();

    enum TableMode {
        Range { start_idx: usize, end_idx: usize },
        Position { pos_idx: usize },
    }

    let mode = match (has_start, has_end, has_pos) {
        (Some(s), Some(e), None) => TableMode::Range {
            start_idx: s,
            end_idx: e,
        },
        (None, None, Some(p)) => TableMode::Position { pos_idx: p },
        (Some(_), Some(_), Some(_)) => {
            return Err(anyhow!(
                "Table has both START/END and POS columns — ambiguous. Use either START+END (range mode) or POS (position mode)."
            ));
        }
        (Some(_), None, _) | (None, Some(_), _) => {
            return Err(anyhow!(
                "Table has only one of START/END — both are required for range mode (found: {:?})",
                headers
            ));
        }
        _ => {
            return Err(anyhow!(
                "Table must have either START+END or POS columns (found: {:?})",
                headers
            ));
        }
    };

    // Validate flank usage
    if let TableMode::Position { .. } = &mode {
        if flank.is_none() {
            return Err(anyhow!(
                "--flank is required when table uses POS column (position mode)"
            ));
        }
    }

    let flank = flank.unwrap_or(0);
    let mut out = Vec::new();

    for (i, rec) in rdr.records().enumerate() {
        let rec = rec?;
        let row = i + 2; // 1-based, accounting for header

        let chr = rec
            .get(chrom_idx)
            .ok_or_else(|| anyhow!("Missing CHROM at row {}", row))?
            .trim()
            .to_string();

        let name = name_idx.and_then(|idx| {
            rec.get(idx)
                .map(|v| v.trim().to_string())
                .filter(|v| !v.is_empty())
        });

        let region = match &mode {
            TableMode::Range {
                start_idx,
                end_idx,
            } => {
                let start = parse_u64_field(&rec, *start_idx, "START", row)?;
                let end = parse_u64_field(&rec, *end_idx, "END", row)?;
                Region {
                    name,
                    chr,
                    start: min(start, end),
                    end: max(start, end),
                }
            }
            TableMode::Position { pos_idx } => {
                let pos = parse_u64_field(&rec, *pos_idx, "POS", row)?;
                Region {
                    name,
                    chr,
                    start: pos.saturating_sub(flank),
                    end: pos + flank,
                }
            }
        };

        out.push(region);
    }

    Ok(out)
}

/// Parse SV CSV/TSV table using column names:
/// CHROM_LEFT, START_LEFT, END_LEFT, POS_LEFT, CHROM_RIGHT, START_RIGHT, END_RIGHT, POS_RIGHT, NAME
fn parse_regions_sv_table(path: &Path, flank: Option<u64>) -> Result<Vec<Region>> {
    let delim = detect_delimiter(path);
    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(true)
        .delimiter(delim)
        .from_path(path)
        .with_context(|| format!("Cannot open SV table: {}", path.display()))?;

    let headers = rdr.headers()?.clone();
    let hmap = build_header_map(&headers);

    // Required chromosome columns
    let chrom_left_idx = hmap.get("CHROM_LEFT").copied().ok_or_else(|| {
        anyhow!(
            "SV table missing required CHROM_LEFT column (found: {:?})",
            headers
        )
    })?;
    let chrom_right_idx = hmap.get("CHROM_RIGHT").copied().ok_or_else(|| {
        anyhow!(
            "SV table missing required CHROM_RIGHT column (found: {:?})",
            headers
        )
    })?;

    // Detect mode from available columns
    let has_start_left = hmap.get("START_LEFT").copied();
    let has_end_left = hmap.get("END_LEFT").copied();
    let has_start_right = hmap.get("START_RIGHT").copied();
    let has_end_right = hmap.get("END_RIGHT").copied();
    let has_pos_left = hmap.get("POS_LEFT").copied();
    let has_pos_right = hmap.get("POS_RIGHT").copied();
    let name_idx = hmap.get("NAME").copied();

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

    let mode = match (
        has_start_left,
        has_end_left,
        has_start_right,
        has_end_right,
        has_pos_left,
        has_pos_right,
    ) {
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
                "SV table must have either START_LEFT+END_LEFT+START_RIGHT+END_RIGHT (range mode) \
                 or POS_LEFT+POS_RIGHT (position mode). Found: {:?}",
                headers
            ));
        }
    };

    if let SvMode::Position { .. } = &mode {
        if flank.is_none() {
            return Err(anyhow!(
                "--flank is required when SV table uses POS_LEFT/POS_RIGHT columns (position mode)"
            ));
        }
    }

    let flank = flank.unwrap_or(0);
    let mut out = Vec::new();

    for (i, rec) in rdr.records().enumerate() {
        let rec = rec?;
        let row = i + 2;

        let chr_left = rec
            .get(chrom_left_idx)
            .ok_or_else(|| anyhow!("Missing CHROM_LEFT at row {}", row))?
            .trim()
            .to_string();
        let chr_right = rec
            .get(chrom_right_idx)
            .ok_or_else(|| anyhow!("Missing CHROM_RIGHT at row {}", row))?
            .trim()
            .to_string();

        let name = name_idx.and_then(|idx| {
            rec.get(idx)
                .map(|v| v.trim().to_string())
                .filter(|v| !v.is_empty())
        });

        match &mode {
            SvMode::Range {
                start_left_idx,
                end_left_idx,
                start_right_idx,
                end_right_idx,
            } => {
                let sl = parse_u64_field(&rec, *start_left_idx, "START_LEFT", row)?;
                let el = parse_u64_field(&rec, *end_left_idx, "END_LEFT", row)?;
                let sr = parse_u64_field(&rec, *start_right_idx, "START_RIGHT", row)?;
                let er = parse_u64_field(&rec, *end_right_idx, "END_RIGHT", row)?;
                out.push(Region {
                    name: name.clone(),
                    chr: chr_left,
                    start: min(sl, el),
                    end: max(sl, el),
                });
                out.push(Region {
                    name,
                    chr: chr_right,
                    start: min(sr, er),
                    end: max(sr, er),
                });
            }
            SvMode::Position {
                pos_left_idx,
                pos_right_idx,
            } => {
                let pl = parse_u64_field(&rec, *pos_left_idx, "POS_LEFT", row)?;
                let pr = parse_u64_field(&rec, *pos_right_idx, "POS_RIGHT", row)?;
                out.push(Region {
                    name: name.clone(),
                    chr: chr_left,
                    start: pl.saturating_sub(flank),
                    end: pl + flank,
                });
                out.push(Region {
                    name,
                    chr: chr_right,
                    start: pr.saturating_sub(flank),
                    end: pr + flank,
                });
            }
        }
    }

    Ok(out)
}

/// Wrap sequence for FASTA output
fn wrap_fasta(seq: &str, width: usize) -> String {
    if seq.is_empty() {
        return String::new();
    }
    let mut out = String::with_capacity(seq.len() + seq.len() / width + 8);
    let mut i = 0usize;
    while i < seq.len() {
        let end = min(i + width, seq.len());
        out.push_str(&seq[i..end]);
        out.push('\n');
        i = end;
    }
    out.trim_end_matches('\n').to_string()
}

fn unique_join(v: &[String], sep: &str) -> String {
    use std::collections::BTreeSet;
    let set: BTreeSet<_> = v.iter().cloned().collect();
    set.into_iter().collect::<Vec<_>>().join(sep)
}

/// Write sequences either combined or per-file
fn write_sequences(
    fasta_path: &Path,
    fai_index: &HashMap<String, FaiRecord>,
    regions: &[Region],
    output_file: Option<&Path>,
    output_dir: Option<&Path>,
    sv_table: Option<&Path>,
) -> Result<()> {
    let sequences: Vec<(Region, String)> = regions
        .par_iter()
        .map(|r| {
            let seq = extract_region(fasta_path, fai_index, r)?;
            let header = format!(">{}:{}-{}\n", r.chr, r.start, r.end);
            Ok::<_, anyhow::Error>((r.clone(), format!("{}{}\n", header, wrap_fasta(&seq, 60))))
        })
        .collect::<Result<_>>()?;

    if let Some(dir) = output_dir {
        std::fs::create_dir_all(dir)?;

        if sv_table.is_some() {
            // Check if number of sequences is even
            if sequences.len() % 2 != 0 {
                return Err(anyhow!(
                    "Number of sequences ({}) is odd; SV table extraction requires an even number of regions (pairs).",
                    sequences.len()
                ));
            }
            // Group sequences into pairs (2 per file)
            let chunked: Vec<_> = sequences.chunks(2).collect();
            chunked
                .into_par_iter()
                .try_for_each(|pair| -> Result<()> {
                    // Use names if present, else chr/start/end for both
                    let filename = if pair.len() == 2 {
                        let (r1, _) = &pair[0];
                        let (r2, _) = &pair[1];
                        match (&r1.name, &r2.name) {
                            (Some(n1), Some(n2)) => {
                                if n1 != n2 {
                                    return Err(anyhow!("Mismatched names in SV pair: '{}' vs '{}'", n1, n2));
                                }
                                format!("{}_{}_{}_{}_{}_{}_{}.fa", n1, r1.chr, r1.start, r1.end, r2.chr, r2.start, r2.end)
                            }
                            _ => format!(
                                "{}_{}_{}_{}_{}_{}.fa",
                                r1.chr, r1.start, r1.end, r2.chr, r2.start, r2.end
                            ),
                        }
                    } else {
                        // Last file if odd count
                        let (r, _) = &pair[0];
                        match &r.name {
                            Some(n) => format!("{}.fa", n),
                            None => format!("{}_{}_{}.fa", r.chr, r.start, r.end),
                        }
                    };
                    let filepath = dir.join(filename);
                    let mut f = File::create(filepath)?;
                    for (_, chunk) in &*pair {
                        f.write_all(chunk.as_bytes())?;
                    }
                    Ok(())
                })?;
        } else {
            sequences
                .into_par_iter()
                .try_for_each(|(region, chunk)| -> Result<()> {
                    let filename = format!("{}_{}_{}.fa", region.chr, region.start, region.end);
                    let filepath = dir.join(filename);
                    let mut f = File::create(filepath)?;
                    f.write_all(chunk.as_bytes())?;
                    Ok(())
                })?;
        }
    } else {
        let mut writer: Box<dyn Write> = match output_file {
            Some(path) => Box::new(BufWriter::new(
                File::options()
                    .create(true)
                    .write(true)
                    .truncate(true)
                    .open(path)?,
            )),
            None => Box::new(BufWriter::new(io::stdout())),
        };
        for (_region, chunk) in sequences {
            writer.write_all(chunk.as_bytes())?;
        }
        writer.flush()?;
    }

    Ok(())
}
