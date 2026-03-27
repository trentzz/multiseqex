use anyhow::{Result, anyhow};
use clap::Parser;
use std::path::PathBuf;

use multiseqex::fai::{build_fai, check_fai_staleness, fai_path_for, read_fai};
use multiseqex::gff::parse_regions_gff;
use multiseqex::intervals::{intersect_regions, subtract_regions};
use multiseqex::mask::MaskIndex;
use multiseqex::output::{OutputConfig, write_sequences};
use multiseqex::region::{
    deduplicate_regions, merge_regions, parse_regions_bed, parse_regions_inline,
    parse_regions_list, resolve_flanks, sort_regions, tile_regions,
};
use multiseqex::stats::write_stats;
use multiseqex::table::{parse_regions_sv_table, parse_regions_table};
use multiseqex::transform::TransformConfig;
use multiseqex::validate::{resolve_bgzip, validate_and_clamp_regions};
use multiseqex::vcf::parse_regions_vcf;

#[derive(Parser, Debug)]
#[command(
    name = "multiseqex",
    author,
    version,
    about = "Multi-sequence extractor for FASTA using FAI"
)]
struct Cli {
    /// Reference FASTA file(s). Supports plain text and bgzip/gzip compressed
    /// files (decompressed transparently). Use `-` for stdin (requires --no-index).
    /// When multiple files are given, contigs are looked up across all files.
    /// A contig must appear in exactly one file.
    #[arg(required = true)]
    fasta: Vec<PathBuf>,

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

    /// Suppress progress messages, warnings, and the progress bar on stderr.
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

    // ── Masking flags ──────────────────────────────────────────────────
    /// BED file defining regions to mask within extracted sequences.
    #[arg(long)]
    mask_bed: Option<PathBuf>,

    /// Replace masked bases with N (default when --mask-bed is given).
    #[arg(long, requires = "mask_bed", conflicts_with = "soft_mask")]
    hard_mask: bool,

    /// Lowercase masked bases instead of replacing with N.
    #[arg(long, requires = "mask_bed", conflicts_with = "hard_mask")]
    soft_mask: bool,

    // ── Transform flags ────────────────────────────────────────────────
    /// Convert T to U in output sequences (DNA to RNA).
    #[arg(long, conflicts_with = "translate")]
    to_rna: bool,

    /// Translate to amino acids (standard genetic code, reading frame from
    /// position 1). Stop codons are represented as *.
    #[arg(long, conflicts_with = "to_rna")]
    translate: bool,

    /// Force all output bases to uppercase.
    #[arg(long, conflicts_with = "lowercase")]
    uppercase: bool,

    /// Force all output bases to lowercase.
    #[arg(long, conflicts_with = "uppercase")]
    lowercase: bool,

    // ── No-index / streaming mode ─────────────────────────────────────
    /// Scan the FASTA sequentially without requiring an FAI index. Loads
    /// the entire FASTA into memory. Allows `-` as the FASTA path for
    /// stdin. Conflicts with --no-build-fai.
    #[arg(long, conflicts_with = "no_build_fai")]
    no_index: bool,

    // ── Interval operations ───────────────────────────────────────────
    /// BED file of intervals to subtract from input regions. Removes
    /// overlapping portions, potentially splitting regions.
    #[arg(long)]
    subtract: Option<PathBuf>,

    /// BED file of intervals to intersect with input regions. Keeps only
    /// overlapping portions.
    #[arg(long)]
    intersect: Option<PathBuf>,

    // ── K-mer tiling ──────────────────────────────────────────────────
    /// Tile each region into windows of this size (in bases).
    #[arg(long)]
    tile: Option<u64>,

    /// Step size for tiling (default: same as --tile, i.e. non-overlapping).
    /// Requires --tile.
    #[arg(long, requires = "tile")]
    step: Option<u64>,

    // ── Alt-seq flags ─────────────────────────────────────────────────
    /// Generate alternate-allele sequences. For each VCF record or table
    /// row, extract the reference context and replace the REF allele with
    /// the ALT allele. Multi-allelic sites produce one output per ALT
    /// allele. Requires --vcf or --table with REF/ALT columns.
    #[arg(long, conflicts_with = "sv_table")]
    alt_seq: bool,

    /// Output both reference and alternate sequences for each variant.
    /// Implies --alt-seq.
    #[arg(long, conflicts_with = "sv_table")]
    alt_seq_both: bool,
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

    // --to-rna and --translate conflict with --stats.
    if (cli.to_rna || cli.translate) && cli.stats {
        return Err(anyhow!(
            "--to-rna and --translate cannot be used with --stats"
        ));
    }

    // --tile must be > 0.
    if let Some(tile) = cli.tile {
        if tile == 0 {
            return Err(anyhow!("--tile must be > 0"));
        }
    }

    // --step must be > 0.
    if let Some(step) = cli.step {
        if step == 0 {
            return Err(anyhow!("--step must be > 0"));
        }
    }

    // --tile conflicts with --sv-table --output-dir (breaks pairing).
    if cli.tile.is_some() && cli.sv_table.is_some() && cli.output_dir.is_some() {
        return Err(anyhow!(
            "--tile cannot be used with --sv-table --output-dir because it breaks \
             the paired region invariant"
        ));
    }

    // --alt-seq-both implies --alt-seq.
    let alt_seq = cli.alt_seq || cli.alt_seq_both;
    let alt_seq_both = cli.alt_seq_both;

    // --alt-seq requires --vcf or --table.
    if alt_seq && cli.vcf.is_none() && cli.table.is_none() {
        return Err(anyhow!(
            "--alt-seq requires --vcf or --table (with REF/ALT columns)"
        ));
    }

    // --no-index requires exactly one FASTA file.
    if cli.no_index && cli.fasta.len() > 1 {
        return Err(anyhow!(
            "--no-index supports only a single FASTA file (or stdin via '-')"
        ));
    }

    // Build mask index if --mask-bed is given.
    let mask_index = if let Some(ref mask_path) = cli.mask_bed {
        Some(MaskIndex::from_bed(mask_path)?)
    } else {
        None
    };

    let mask_mode = if cli.soft_mask {
        multiseqex::mask::MaskMode::Soft
    } else {
        // Default to hard mask when --mask-bed is given.
        multiseqex::mask::MaskMode::Hard
    };

    // Build transform config.
    let transform = TransformConfig {
        to_rna: cli.to_rna,
        uppercase: cli.uppercase,
        lowercase: cli.lowercase,
        translate: cli.translate,
    };

    // ── Load FASTA files and build combined index ──────────────────────

    // Validate all FASTA files and build/load their FAI indices.
    // Bgzip/gzip files are decompressed to temporary files transparently.
    let mut fasta_paths: Vec<PathBuf> = Vec::with_capacity(cli.fasta.len());
    let mut fai_indices = Vec::new();
    // Keep temp directories alive for the duration of the programme.
    let mut _bgzip_temps: Vec<tempfile::TempDir> = Vec::new();
    // In-memory sequences for --no-index mode.
    let mut noindex_sequences: Option<std::collections::HashMap<String, Vec<u8>>> = None;

    let is_stdin = cli.fasta.len() == 1 && cli.fasta[0] == std::path::Path::new("-");

    if cli.no_index {
        // --no-index: load the FASTA into memory without an FAI index.
        if !cli.quiet {
            if is_stdin {
                eprintln!("Reading FASTA from stdin (--no-index mode, loading into memory).");
            } else {
                eprintln!(
                    "Loading FASTA into memory (--no-index mode): {}",
                    cli.fasta[0].display()
                );
            }
        }
        let reader: Box<dyn std::io::Read> = if is_stdin {
            Box::new(std::io::stdin())
        } else {
            let fasta_path = &cli.fasta[0];
            if multiseqex::validate::is_gzip(fasta_path)? {
                let f = std::fs::File::open(fasta_path)?;
                Box::new(flate2::read::GzDecoder::new(f))
            } else {
                Box::new(std::fs::File::open(fasta_path)?)
            }
        };
        let (seqs, fai) = multiseqex::noindex::load_fasta_into_memory(reader)?;
        // We still need a dummy fasta_paths entry for the output pipeline.
        // In no-index mode, we write output directly and do not use the file-based
        // extraction pipeline, so a placeholder suffices.
        fasta_paths.push(if is_stdin {
            PathBuf::from("-")
        } else {
            cli.fasta[0].clone()
        });
        fai_indices.push(fai);
        noindex_sequences = Some(seqs);
    } else {
        if is_stdin {
            return Err(anyhow!("Reading from stdin ('-') requires --no-index"));
        }
        for fasta_path in &cli.fasta {
            let (resolved_path, tmp_dir) = resolve_bgzip(fasta_path, cli.quiet)?;
            if let Some(td) = tmp_dir {
                _bgzip_temps.push(td);
            }
            let fai_path = fai_path_for(fasta_path);
            // For bgzip files, try the FAI next to the original file first.
            // If not found, also try next to the decompressed temp file.
            let effective_fai = if fai_path.exists() {
                fai_path.clone()
            } else {
                fai_path_for(&resolved_path)
            };
            let mut fai_just_built = false;
            if !effective_fai.exists() {
                if cli.no_build_fai {
                    return Err(anyhow!(
                        "Missing index: {} (use samtools faidx or remove --no-build-fai)",
                        fai_path.display()
                    ));
                }
                if !cli.quiet {
                    eprintln!("Index not found. Building FAI: {}", effective_fai.display());
                }
                build_fai(&resolved_path, &effective_fai)?;
                fai_just_built = true;
            }

            if !fai_just_built && !cli.quiet {
                check_fai_staleness(&resolved_path, &effective_fai);
            }
            let idx = read_fai(&effective_fai)?;
            fasta_paths.push(resolved_path);
            fai_indices.push(idx);
        }
    }

    // Build the combined FAI index. Maps contig name to (fasta_file_index, FaiRecord).
    // Error if a contig appears in multiple FASTA files.
    let mut combined_fai = std::collections::HashMap::new();
    let mut contig_to_fasta: std::collections::HashMap<String, usize> =
        std::collections::HashMap::new();
    for (file_idx, idx) in fai_indices.iter().enumerate() {
        for (name, rec) in idx {
            if let Some(prev_idx) = contig_to_fasta.get(name) {
                return Err(anyhow!(
                    "Contig '{}' appears in multiple FASTA files: '{}' and '{}'",
                    name,
                    fasta_paths[*prev_idx].display(),
                    fasta_paths[file_idx].display()
                ));
            }
            contig_to_fasta.insert(name.clone(), file_idx);
            combined_fai.insert(name.clone(), rec.clone());
        }
    }

    let fai_index = combined_fai;

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
        if alt_seq {
            // Validate and use the alt-seq table parser.
            multiseqex::table::validate_table_alt_seq(p, cli.delimiter.as_deref())?;
            let base_idx = regions.len();
            let table_regions = multiseqex::table::parse_regions_table_alt_seq(
                p,
                cli.flank,
                cli.flank_left,
                cli.flank_right,
                cli.delimiter.as_deref(),
                alt_seq_both,
            )?;
            for r in &table_regions {
                // Pad vcf_descriptions to align with region indices.
                while vcf_descriptions.len() < base_idx {
                    vcf_descriptions.push(String::new());
                }
                if let Some(ref info) = r.alt_info {
                    let mut desc = format!("REF={} ALT={}", info.ref_allele, info.alt_allele);
                    if alt_seq_both {
                        desc.push_str(" alt_seq");
                    }
                    vcf_descriptions.push(desc);
                } else if alt_seq_both {
                    // Reference entry in both mode. We don't know the ALT yet
                    // but next entry will have it. Use a placeholder from the
                    // next region's alt_info if available.
                    vcf_descriptions.push("ref_seq".to_string());
                } else {
                    vcf_descriptions.push(String::new());
                }
            }
            regions.extend(table_regions);
        } else {
            regions.extend(parse_regions_table(
                p,
                cli.flank,
                cli.flank_left,
                cli.flank_right,
                cli.delimiter.as_deref(),
            )?);
        }
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

        // When --alt-seq is active, expand multi-allelic sites and attach
        // AltInfo to each region so the output pipeline performs the
        // substitution.
        let vcf_results = if alt_seq {
            multiseqex::vcf::expand_vcf_alt_seq(vcf_results, alt_seq_both)
        } else {
            vcf_results
        };

        let base_idx = regions.len();
        for (region, rec) in &vcf_results {
            regions.push(region.clone());
            // Pad vcf_descriptions to align with region indices.
            while vcf_descriptions.len() < base_idx {
                vcf_descriptions.push(String::new());
            }
            let mut desc = multiseqex::vcf::vcf_description(rec);
            if alt_seq_both {
                if region.alt_info.is_some() {
                    desc.push_str(" alt_seq");
                } else {
                    desc.push_str(" ref_seq");
                }
            }
            vcf_descriptions.push(desc);
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
                alt_info: None,
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
                alt_info: None,
            });
        }
    }

    if regions.is_empty() {
        return Err(anyhow!(
            "No regions provided. Use --regions, --list, --bed, --table, --sv-table, \
             --vcf, --gff, --contigs, or --contig-list."
        ));
    }

    // ── Interval operations (before dedup/sort/merge) ────────────────
    // Intersect first, then subtract (as specified).
    if let Some(ref bed_path) = cli.intersect {
        regions = intersect_regions(&regions, bed_path)?;
        if regions.is_empty() {
            return Err(anyhow!(
                "No regions remain after --intersect. All input regions were outside \
                 the intervals in '{}'.",
                bed_path.display()
            ));
        }
        if !cli.quiet {
            eprintln!("After --intersect: {} region(s) remain.", regions.len());
        }
    }
    if let Some(ref bed_path) = cli.subtract {
        regions = subtract_regions(&regions, bed_path)?;
        if regions.is_empty() {
            return Err(anyhow!(
                "No regions remain after --subtract. All input regions were fully \
                 covered by intervals in '{}'.",
                bed_path.display()
            ));
        }
        if !cli.quiet {
            eprintln!("After --subtract: {} region(s) remain.", regions.len());
        }
    }

    // ── K-mer tiling ─────────────────────────────────────────────────
    if let Some(tile_size) = cli.tile {
        let step = cli.step.unwrap_or(tile_size);
        let before = regions.len();
        regions = tile_regions(&regions, tile_size, step);
        if !cli.quiet {
            eprintln!(
                "Tiled {before} region(s) into {} tiles (tile={tile_size}, step={step}).",
                regions.len()
            );
        }
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
        if noindex_sequences.is_some() {
            return Err(anyhow!("--stats is not supported with --no-index"));
        }
        write_stats(&fasta_paths, &fai_index, &contig_to_fasta, &regions)?;
        return Ok(());
    }

    let effective_line_width = if cli.no_wrap { 0 } else { cli.line_width };

    let config = OutputConfig {
        fastq: cli.fastq,
        qual_char: cli.qual.chars().next().unwrap_or('I'),
        name_template: cli.name_template,
        vcf_descriptions,
        mask_index,
        mask_mode,
        transform,
    };

    let is_sv = cli.sv_table.is_some();
    let output_to_file = cli.output.is_some() || cli.output_dir.is_some();

    // Show progress bar on stderr when writing to a file and not quiet.
    let show_progress = !cli.quiet && output_to_file;

    // ── --no-index extraction path ───────────────────────────────────
    if let Some(ref sequences) = noindex_sequences {
        use multiseqex::extract::reverse_complement;
        use multiseqex::noindex::extract_region_from_memory;
        use multiseqex::output::wrap_fasta;

        if show_progress {
            eprintln!("Extracting {} regions (--no-index mode)...", regions.len());
        }
        let mut writer: Box<dyn std::io::Write> = match cli.output.as_deref() {
            Some(p) => Box::new(std::io::BufWriter::new(std::fs::File::create(p)?)),
            None => Box::new(std::io::BufWriter::new(std::io::stdout())),
        };
        for (i, r) in regions.iter().enumerate() {
            let mut seq = extract_region_from_memory(sequences, r)?;
            // Apply alt-seq substitution if applicable.
            if let Some(ref info) = r.alt_info {
                seq = multiseqex::extract::apply_alt_substitution(&seq, r.start, info)?;
            }
            // Apply strand and RC.
            let strand_rc = r.strand == Some('-');
            if strand_rc ^ cli.reverse_complement {
                seq = reverse_complement(&seq);
            }
            // Apply masking and transforms.
            if let Some(ref mask_idx) = config.mask_index {
                seq = mask_idx.apply(&seq, &r.chr, r.start, r.end, config.mask_mode);
            }
            if config.transform.any_active() {
                seq = config.transform.apply(&seq);
            }
            // Format output.
            if cli.tab_out {
                let name = r.name.as_deref().unwrap_or(".");
                writeln!(
                    writer,
                    "{}\t{}\t{}\t{}\t{}",
                    r.chr, r.start, r.end, name, seq
                )?;
            } else if config.fastq {
                let qual: String = std::iter::repeat_n(config.qual_char, seq.len()).collect();
                let header = format_header(r, i, &config);
                write!(writer, "@{header}\n{seq}\n+\n{qual}\n")?;
            } else {
                let header = format_header(r, i, &config);
                write!(
                    writer,
                    ">{header}\n{}\n",
                    wrap_fasta(&seq, effective_line_width)
                )?;
            }
        }
        use std::io::Write as _;
        writer.flush()?;
        if show_progress {
            eprintln!("Done.");
        }
        return Ok(());
    }

    // ── Standard FAI-based extraction path ───────────────────────────
    if show_progress {
        eprintln!("Extracting {} regions from FASTA...", regions.len());
    }
    write_sequences(
        &fasta_paths,
        &fai_index,
        &contig_to_fasta,
        &regions,
        cli.output.as_deref(),
        cli.output_dir.as_deref(),
        is_sv,
        cli.reverse_complement,
        effective_line_width,
        cli.tab_out,
        &config,
        show_progress,
    )?;
    if show_progress {
        eprintln!("Done.");
    }
    Ok(())
}

/// Format a FASTA/FASTQ header for a region.
fn format_header(r: &multiseqex::region::Region, index: usize, config: &OutputConfig) -> String {
    if let Some(ref template) = config.name_template {
        multiseqex::template::expand_template(template, r, index + 1)
    } else {
        let suffix = match r.strand {
            Some('+') => "(+)",
            Some('-') => "(-)",
            Some('.') => "(.)",
            _ => "",
        };
        match &r.name {
            Some(name) => format!("{name} {}:{}-{}{suffix}", r.chr, r.start, r.end),
            None => format!("{}:{}-{}{suffix}", r.chr, r.start, r.end),
        }
    }
}
