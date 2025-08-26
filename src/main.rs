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

    /// CSV/TSV with 3 columns: chr,start,end (header allowed)
    /// Delimiter auto-detected from extension: .tsv => tab, else comma
    #[arg(long)]
    table: Option<PathBuf>,

    /// Flank size (only used if table has chr,pos with 2 columns)
    #[arg(long, default_value_t = 0)]
    flank: u64,

    /// Output FASTA file (single combined; default: stdout)
    #[arg(short, long)]
    output: Option<PathBuf>,

    /// Output directory for per-sequence FASTA files
    #[arg(long)]
    output_dir: Option<PathBuf>,

    /// Number of worker threads (default: Rayon default, usually #cpus)
    #[arg(long)]
    threads: Option<usize>,

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
        regions.extend(parse_regions_list(p,  cli.flank)?);
    }
    if let Some(p) = cli.table.as_ref() {
        regions.extend(parse_regions_table(p, cli.flank)?);
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
fn parse_regions_inline(s: &str, flank: u64) -> Result<Vec<Region>> {
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
fn parse_regions_list(path: &Path, flank: u64) -> Result<Vec<Region>> {
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
fn parse_region_str(s: &str, flank: u64) -> Result<Region> {
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

    let pos: u64 = pos.replace(',', "").parse()?;
    let flank_val: u64 = flank_val.replace(',', "").parse()?;
    let start = pos.saturating_sub(flank_val.max(flank));
    let end = pos + flank_val.max(flank);

    Ok(Region {
        chr: chr.to_string(),
        start: start.max(1),
        end,
    })
}

/// Parse CSV/TSV table (2 or 3 columns)
fn parse_regions_table(path: &Path, flank: u64) -> Result<Vec<Region>> {
    let delim = if path
        .extension()
        .and_then(|e| e.to_str())
        .map(|e| e.eq_ignore_ascii_case("tsv"))
        .unwrap_or(false)
    {
        b'\t'
    } else {
        b','
    };
    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(true)
        .delimiter(delim)
        .from_path(path)
        .with_context(|| format!("Cannot open table: {}", path.display()))?;
    let headers = rdr.headers().map(|h| h.clone()).ok();

    let mut out = Vec::new();
    let mut records = rdr.records();

    // Determine mode from first record
    let first = records
        .next()
        .ok_or_else(|| anyhow!("Empty table: {}", path.display()))??;
    let mode = first.len();
    if mode != 2 && mode != 3 {
        return Err(anyhow!("Table must have 2 or 3 columns (got {})", mode));
    }

    let process_record = |rec: csv::StringRecord, idx: usize| -> Result<Region> {
        match mode {
            2 => {
                let chr = rec.get(0).unwrap().trim();
                let pos: u64 = rec
                    .get(1)
                    .unwrap()
                    .trim()
                    .replace(',', "")
                    .parse()
                    .with_context(|| {
                        format!("Bad pos at row {} (headers {:?})", idx + 1, headers)
                    })?;
                let start = pos.saturating_sub(flank);
                let end = pos + flank;
                Ok(Region {
                    chr: chr.to_string(),
                    start,
                    end,
                })
            }
            3 => {
                let chr = rec.get(0).unwrap().trim();
                let start: u64 = rec
                    .get(1)
                    .unwrap()
                    .trim()
                    .replace(',', "")
                    .parse()
                    .with_context(|| {
                        format!("Bad start at row {} (headers {:?})", idx + 1, headers)
                    })?;
                let end: u64 = rec
                    .get(2)
                    .unwrap()
                    .trim()
                    .replace(',', "")
                    .parse()
                    .with_context(|| {
                        format!("Bad end at row {} (headers {:?})", idx + 1, headers)
                    })?;
                Ok(Region {
                    chr: chr.to_string(),
                    start: min(start, end),
                    end: max(start, end),
                })
            }
            _ => unreachable!(),
        }
    };

    out.push(process_record(first, 0)?);
    for (i, rec) in records.enumerate() {
        let rec = rec?;
        if rec.len() != mode {
            return Err(anyhow!(
                "Row {} has {} columns but expected {}",
                i + 2,
                rec.len(),
                mode
            ));
        }
        out.push(process_record(rec, i + 1)?);
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
        sequences
            .into_par_iter()
            .try_for_each(|(region, chunk)| -> Result<()> {
                let filename = format!("{}_{}_{}.fa", region.chr, region.start, region.end);
                let filepath = dir.join(filename);
                let mut f = File::create(filepath)?;
                f.write_all(chunk.as_bytes())?;
                Ok(())
            })?;
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
