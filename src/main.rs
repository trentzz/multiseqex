mod extract;
mod fai;
mod output;
mod region;
mod table;
mod validate;

use anyhow::{Result, anyhow};
use clap::Parser;
use std::path::PathBuf;

use fai::{build_fai, check_fai_staleness, fai_path_for, read_fai};
use output::write_sequences;
use region::{parse_regions_inline, parse_regions_list};
use table::{parse_regions_sv_table, parse_regions_table};
use validate::{detect_gzip_and_reject, validate_and_clamp_regions};

#[derive(Parser, Debug)]
#[command(
    name = "multiseqex",
    author,
    version,
    about = "Multi-sequence extractor for FASTA using FAI"
)]
struct Cli {
    /// Reference FASTA file (plain text, not compressed).
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
    /// Extra columns are ignored. Delimiter auto-detected (.tsv → tab, else sniff
    /// for tabs, else comma). Use --delimiter to override.
    #[arg(long)]
    table: Option<PathBuf>,

    /// CSV/TSV structural-variant table with named columns.
    ///
    /// Required: CHROM_LEFT, CHROM_RIGHT.
    /// Range mode: START_LEFT, END_LEFT, START_RIGHT, END_RIGHT.
    /// Position mode: POS_LEFT, POS_RIGHT (requires --flank).
    /// Optional: NAME.
    /// Extra columns are ignored.
    #[arg(long, conflicts_with_all = ["regions", "table", "list"])]
    sv_table: Option<PathBuf>,

    /// Flank size for position-mode tables (required with POS columns).
    #[arg(long)]
    flank: Option<u64>,

    /// Output FASTA file (single combined file; default: stdout).
    #[arg(short, long, conflicts_with = "output_dir")]
    output: Option<PathBuf>,

    /// Output directory, one FASTA file per region (or per SV pair).
    #[arg(long, conflicts_with = "output")]
    output_dir: Option<PathBuf>,

    /// Number of worker threads (default: all available CPUs).
    #[arg(long)]
    threads: Option<usize>,

    /// Override delimiter for --table / --sv-table.
    ///
    /// Accepts "tab", "comma", or a single character (e.g. ";").
    /// When omitted the delimiter is auto-detected: .tsv → tab,
    /// otherwise the first line is sniffed for tabs, falling back to comma.
    #[arg(long)]
    delimiter: Option<String>,

    /// Error if .fai is missing instead of building one automatically.
    #[arg(long)]
    no_build_fai: bool,
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    if let Some(t) = cli.threads
        && let Err(e) = rayon::ThreadPoolBuilder::new()
            .num_threads(t)
            .build_global()
    {
        eprintln!("Warning: failed to build thread pool with {t} threads: {e}");
    }

    detect_gzip_and_reject(&cli.fasta)?;
    let fai_path = fai_path_for(&cli.fasta);
    let mut fai_just_built = false;
    if !fai_path.exists() {
        if cli.no_build_fai {
            return Err(anyhow!(
                "Missing index: {} (use samtools faidx or remove --no-build-fai)",
                fai_path.display()
            ));
        }
        eprintln!("Index not found. Building FAI: {}", fai_path.display());
        build_fai(&cli.fasta, &fai_path)?;
        fai_just_built = true;
    }

    if !fai_just_built {
        check_fai_staleness(&cli.fasta, &fai_path);
    }
    let fai_index = read_fai(&fai_path)?;
    let mut regions = Vec::new();
    if let Some(s) = cli.regions.as_deref() {
        regions.extend(parse_regions_inline(s, cli.flank)?);
    }
    if let Some(p) = cli.list.as_ref() {
        regions.extend(parse_regions_list(p, cli.flank)?);
    }
    if let Some(p) = cli.table.as_ref() {
        regions.extend(parse_regions_table(p, cli.flank, cli.delimiter.as_deref())?);
    }
    if let Some(p) = cli.sv_table.as_ref() {
        regions.extend(parse_regions_sv_table(
            p,
            cli.flank,
            cli.delimiter.as_deref(),
        )?);
    }
    if regions.is_empty() {
        return Err(anyhow!(
            "No regions provided. Use --regions, --list, --table, or --sv-table."
        ));
    }
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
