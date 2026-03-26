use anyhow::{Result, anyhow};
use clap::Parser;
use std::path::PathBuf;

use multiseqex::fai::{build_fai, check_fai_staleness, fai_path_for, read_fai};
use multiseqex::output::write_sequences;
use multiseqex::region::{
    deduplicate_regions, parse_regions_bed, parse_regions_inline, parse_regions_list, sort_regions,
};
use multiseqex::table::{parse_regions_sv_table, parse_regions_table};
use multiseqex::validate::{detect_gzip_and_reject, validate_and_clamp_regions};

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

    /// BED file (tab-separated: chr, start, end, optional name).
    ///
    /// BED uses 0-based half-open coordinates. They are converted to
    /// 1-based inclusive internally (start+1, end unchanged).
    #[arg(long, conflicts_with = "sv_table")]
    bed: Option<PathBuf>,

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

    /// Flank size for position-mode regions (POS columns in tables and
    /// pos+flank inline syntax). Has no effect on BED or range-format regions.
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

    /// Reverse complement all extracted sequences.
    #[arg(long = "rc", alias = "reverse-complement")]
    reverse_complement: bool,

    /// Deduplicate regions with identical chr, start, end before extraction.
    #[arg(long = "dedup", alias = "deduplicate")]
    deduplicate: bool,

    /// Sort output regions by genomic coordinate (natural chromosome order, then start).
    #[arg(long)]
    sort: bool,

    /// Suppress progress messages and warnings on stderr.
    #[arg(short, long)]
    quiet: bool,

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
        && !cli.quiet
    {
        eprintln!("Warning: failed to build thread pool with {t} threads: {e}");
    }

    // --delimiter only makes sense with --table or --sv-table.
    if cli.delimiter.is_some() && cli.table.is_none() && cli.sv_table.is_none() {
        return Err(anyhow!("--delimiter requires --table or --sv-table"));
    }

    // --dedup/--sort with --sv-table --output-dir would break paired region
    // ordering. Forbid this combination.
    if (cli.deduplicate || cli.sort) && cli.sv_table.is_some() && cli.output_dir.is_some() {
        let flag = if cli.deduplicate && cli.sort {
            "--dedup and --sort"
        } else if cli.deduplicate {
            "--dedup"
        } else {
            "--sort"
        };
        return Err(anyhow!(
            "{flag} cannot be used with --sv-table --output-dir because they break \
             the paired region invariant (pairs may be deduplicated or reordered)"
        ));
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
        if !cli.quiet {
            eprintln!("Index not found. Building FAI: {}", fai_path.display());
        }
        build_fai(&cli.fasta, &fai_path)?;
        fai_just_built = true;
    }

    if !fai_just_built && !cli.quiet {
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
    if let Some(p) = cli.bed.as_ref() {
        regions.extend(parse_regions_bed(p, cli.flank)?);
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
            "No regions provided. Use --regions, --list, --bed, --table, or --sv-table."
        ));
    }
    if cli.deduplicate {
        let removed = deduplicate_regions(&mut regions);
        if removed > 0 && !cli.quiet {
            eprintln!("Deduplicated: removed {removed} duplicate region(s).");
        }
    }
    if cli.sort {
        sort_regions(&mut regions);
    }
    validate_and_clamp_regions(&mut regions, &fai_index)?;
    let is_sv = cli.sv_table.is_some();
    let output_to_file = cli.output.is_some() || cli.output_dir.is_some();
    if !cli.quiet && output_to_file {
        eprintln!("Extracting {} regions from FASTA...", regions.len());
    }
    write_sequences(
        &cli.fasta,
        &fai_index,
        &regions,
        cli.output.as_deref(),
        cli.output_dir.as_deref(),
        is_sv,
        cli.reverse_complement,
    )?;
    if !cli.quiet && output_to_file {
        eprintln!("Done.");
    }
    Ok(())
}
