use anyhow::{Result, anyhow};
use clap::Parser;
use std::path::PathBuf;

use multiseqex::fai::{build_fai, check_fai_staleness, fai_path_for, read_fai};
use multiseqex::gff::parse_regions_gff;
use multiseqex::output::{OutputConfig, write_sequences};
use multiseqex::region::{
    deduplicate_regions, merge_regions, parse_regions_bed, parse_regions_inline,
    parse_regions_list, resolve_flanks, sort_regions,
};
use multiseqex::stats::write_stats;
use multiseqex::table::{parse_regions_sv_table, parse_regions_table};
use multiseqex::validate::{detect_gzip_and_reject, validate_and_clamp_regions};
use multiseqex::vcf::parse_regions_vcf;

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
    /// Optional: NAME, STRAND.
    /// Extra columns are ignored. Delimiter auto-detected (.tsv -> tab, else sniff
    /// for tabs, else comma). Use --delimiter to override.
    #[arg(long)]
    table: Option<PathBuf>,

    /// CSV/TSV structural-variant table with named columns.
    ///
    /// Required: CHROM_LEFT, CHROM_RIGHT.
    /// Range mode: START_LEFT, END_LEFT, START_RIGHT, END_RIGHT.
    /// Position mode: POS_LEFT, POS_RIGHT (requires --flank).
    /// Optional: NAME, STRAND.
    /// Extra columns are ignored.
    #[arg(long, conflicts_with_all = ["regions", "table", "list", "contigs", "contig_list", "vcf", "gff"])]
    sv_table: Option<PathBuf>,

    /// VCF file. Each record produces a region spanning POS to POS+len(REF)-1.
    /// The ID field is used as the region name (unless "."). REF and ALT are
    /// included in the FASTA header description.
    #[arg(long, conflicts_with = "sv_table")]
    vcf: Option<PathBuf>,

    /// GFF3/GTF annotation file. Extracts regions for features matching
    /// --gff-feature (default: "gene").
    #[arg(long, conflicts_with = "sv_table")]
    gff: Option<PathBuf>,

    /// Feature type to filter when using --gff (default: "gene").
    #[arg(long, default_value = "gene", requires = "gff")]
    gff_feature: String,

    /// Flank size for position-mode regions (POS columns in tables and
    /// pos+flank inline syntax). Has no effect on BED or range-format regions.
    #[arg(long, conflicts_with_all = ["flank_left", "flank_right"])]
    flank: Option<u64>,

    /// Left-side flank size. Overrides --flank on the left side.
    #[arg(long)]
    flank_left: Option<u64>,

    /// Right-side flank size. Overrides --flank on the right side.
    #[arg(long)]
    flank_right: Option<u64>,

    /// Comma-separated contig names to extract in full.
    #[arg(long, conflicts_with = "sv_table")]
    contigs: Option<String>,

    /// File with one contig name per line to extract in full.
    #[arg(long, conflicts_with = "sv_table")]
    contig_list: Option<PathBuf>,

    /// Output FASTA file (single combined file; default: stdout).
    #[arg(short, long, conflicts_with = "output_dir")]
    output: Option<PathBuf>,

    /// Output directory, one FASTA file per region (or per SV pair).
    #[arg(long, conflicts_with_all = ["output", "tab_out"])]
    output_dir: Option<PathBuf>,

    /// Number of worker threads (default: all available CPUs).
    #[arg(long)]
    threads: Option<usize>,

    /// Override delimiter for --table / --sv-table.
    ///
    /// Accepts "tab", "comma", or a single character (e.g. ";").
    /// When omitted the delimiter is auto-detected: .tsv -> tab,
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

    /// FASTA line width (default 60). Set to 0 to disable wrapping.
    #[arg(long, default_value = "60", conflicts_with = "no_wrap")]
    line_width: usize,

    /// Disable FASTA line wrapping (shorthand for --line-width 0).
    #[arg(long, conflicts_with = "line_width")]
    no_wrap: bool,

    /// Emit TSV output instead of FASTA (columns: chr, start, end, name, sequence).
    #[arg(long, conflicts_with_all = ["output_dir", "fastq"])]
    tab_out: bool,

    /// Merge overlapping or book-ended regions on the same chromosome.
    /// Implies --sort.
    #[arg(long)]
    merge: bool,

    /// Maximum distance between regions to merge (default 0). Requires --merge.
    #[arg(long, default_value = "0")]
    merge_distance: u64,

    /// Emit FASTQ output instead of FASTA. Uses a constant quality character
    /// (default 'I', phred 40). See --qual.
    #[arg(long, conflicts_with_all = ["tab_out"])]
    fastq: bool,

    /// Quality character for FASTQ output (default 'I', phred 40).
    /// Must be a single ASCII character.
    #[arg(long, default_value = "I", requires = "fastq")]
    qual: String,

    /// Print per-region statistics (TSV) instead of extracting sequences.
    /// Columns: chr, start, end, name, length, gc_percent, n_count, masked_count.
    #[arg(long, conflicts_with_all = ["output", "output_dir", "fastq", "tab_out", "reverse_complement"])]
    stats: bool,

    /// Template for FASTA/FASTQ header names. Placeholders: {chr}, {start},
    /// {end}, {name}, {length}, {index}, {strand}.
    #[arg(long, conflicts_with_all = ["tab_out", "stats"])]
    name_template: Option<String>,
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

    // --flank-left/--flank-right must be given together.
    if cli.flank_left.is_some() != cli.flank_right.is_some() {
        return Err(anyhow!(
            "--flank-left and --flank-right must be specified together"
        ));
    }

    // --merge-distance requires --merge.
    if cli.merge_distance > 0 && !cli.merge {
        return Err(anyhow!("--merge-distance requires --merge"));
    }

    // --merge conflicts with --sv-table --output-dir.
    if cli.merge && cli.sv_table.is_some() && cli.output_dir.is_some() {
        return Err(anyhow!(
            "--merge cannot be used with --sv-table --output-dir because it breaks \
             the paired region invariant"
        ));
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

    // Validate --qual is a single ASCII character.
    if cli.qual.len() != 1 || !cli.qual.is_ascii() {
        return Err(anyhow!(
            "--qual must be a single ASCII character, got '{}'",
            cli.qual
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
    let mut vcf_descriptions: Vec<String> = Vec::new();

    if let Some(s) = cli.regions.as_deref() {
        regions.extend(parse_regions_inline(s, cli.flank)?);
    }
    if let Some(p) = cli.list.as_ref() {
        regions.extend(parse_regions_list(p, cli.flank)?);
    }
    if let Some(p) = cli.bed.as_ref() {
        regions.extend(parse_regions_bed(
            p,
            cli.flank,
            cli.flank_left,
            cli.flank_right,
        )?);
    }
    if let Some(p) = cli.table.as_ref() {
        regions.extend(parse_regions_table(
            p,
            cli.flank,
            cli.flank_left,
            cli.flank_right,
            cli.delimiter.as_deref(),
        )?);
    }
    if let Some(p) = cli.sv_table.as_ref() {
        regions.extend(parse_regions_sv_table(
            p,
            cli.flank,
            cli.flank_left,
            cli.flank_right,
            cli.delimiter.as_deref(),
        )?);
    }

    // --vcf: extract regions from VCF records.
    if let Some(p) = cli.vcf.as_ref() {
        let vcf_results = parse_regions_vcf(p, cli.flank, cli.flank_left, cli.flank_right)?;
        let base_idx = regions.len();
        for (region, rec) in vcf_results {
            regions.push(region);
            // Pad vcf_descriptions to align with region indices.
            while vcf_descriptions.len() < base_idx {
                vcf_descriptions.push(String::new());
            }
            vcf_descriptions.push(multiseqex::vcf::vcf_description(&rec));
        }
    }

    // --gff: extract regions from GFF3/GTF annotation.
    if let Some(p) = cli.gff.as_ref() {
        let gff_regions = parse_regions_gff(p, &cli.gff_feature)?;
        // Apply flanking to GFF regions if specified.
        let has_flank =
            cli.flank.is_some() || cli.flank_left.is_some() || cli.flank_right.is_some();
        if has_flank {
            let (fl, fr) = resolve_flanks(cli.flank, cli.flank_left, cli.flank_right);
            for mut r in gff_regions {
                r.start = r.start.saturating_sub(fl).max(1);
                r.end = r.end.saturating_add(fr);
                regions.push(r);
            }
        } else {
            regions.extend(gff_regions);
        }
    }

    // --contigs: extract whole contigs by name.
    if let Some(contig_str) = cli.contigs.as_deref() {
        for name in contig_str.split(',').filter(|s| !s.trim().is_empty()) {
            let name = name.trim();
            let rec = fai_index
                .get(name)
                .ok_or_else(|| anyhow!("Contig '{}' not found in FAI index", name))?;
            regions.push(multiseqex::region::Region {
                name: None,
                chr: name.to_string(),
                start: 1,
                end: rec.length,
                strand: None,
            });
        }
    }

    // --contig-list: extract whole contigs from a file (one per line).
    if let Some(p) = cli.contig_list.as_ref() {
        let content = std::fs::read_to_string(p)
            .map_err(|e| anyhow!("Cannot read contig list file '{}': {}", p.display(), e))?;
        for line in content.lines() {
            let name = line.trim();
            if name.is_empty() || name.starts_with('#') {
                continue;
            }
            let rec = fai_index
                .get(name)
                .ok_or_else(|| anyhow!("Contig '{}' not found in FAI index", name))?;
            regions.push(multiseqex::region::Region {
                name: None,
                chr: name.to_string(),
                start: 1,
                end: rec.length,
                strand: None,
            });
        }
    }

    if regions.is_empty() {
        return Err(anyhow!(
            "No regions provided. Use --regions, --list, --bed, --table, --sv-table, \
             --vcf, --gff, --contigs, or --contig-list."
        ));
    }
    if cli.deduplicate {
        let removed = deduplicate_regions(&mut regions);
        if removed > 0 && !cli.quiet {
            eprintln!("Deduplicated: removed {removed} duplicate region(s).");
        }
    }

    // --merge implies --sort.
    if cli.merge {
        let before = regions.len();
        merge_regions(&mut regions, cli.merge_distance);
        let after = regions.len();
        if before != after && !cli.quiet {
            eprintln!(
                "Merged: {before} regions into {after} ({} merged).",
                before - after
            );
        }
    } else if cli.sort {
        sort_regions(&mut regions);
    }

    validate_and_clamp_regions(&mut regions, &fai_index)?;

    // --stats: print statistics table and exit.
    if cli.stats {
        write_stats(&cli.fasta, &fai_index, &regions)?;
        return Ok(());
    }

    let effective_line_width = if cli.no_wrap { 0 } else { cli.line_width };

    let config = OutputConfig {
        fastq: cli.fastq,
        qual_char: cli.qual.chars().next().unwrap_or('I'),
        name_template: cli.name_template,
        vcf_descriptions,
    };

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
        effective_line_width,
        cli.tab_out,
        &config,
    )?;
    if !cli.quiet && output_to_file {
        eprintln!("Done.");
    }
    Ok(())
}
