# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.1] - 2026-03-28

### Added

- **Alternate allele sequences** (`--alt-seq`, `--alt-seq-both`): generate
  sequences with the ALT allele substituted in place of REF. Works with `--vcf`
  and `--table` (requires REF/ALT columns). Multi-allelic sites produce one
  output per ALT allele. `--alt-seq-both` outputs both the reference and
  alternate sequences for each variant.

### Fixed

- Removed security audit CI job that was failing with "Resource not accessible
  by integration".

## [0.2.0] - 2026-03-27

### Added

- **VCF input** (`--vcf`): extracts regions spanning each VCF record's REF
  allele. The ID field is used as the region name; REF and ALT are included in
  the FASTA header description. Supports flanking via `--flank`,
  `--flank-left`, `--flank-right`.
- **GFF3/GTF input** (`--gff`, `--gff-feature`): extracts regions for features
  matching a given type (default: `gene`). Supports flanking.
- **Whole contig extraction** (`--contigs`, `--contig-list`): extract entire
  contigs by name, either inline (comma-separated) or from a file.
- **BED file input** (`--bed`): 0-based half-open to 1-based inclusive
  coordinate conversion. Supports optional name column (col 4) and strand
  column (col 6). Rejects empty intervals and malformed coordinates.
- **BED flanking**: `--flank` extends BED regions symmetrically after
  coordinate conversion.
- **Asymmetric flanking** (`--flank-left`, `--flank-right`): apply different
  flank sizes to each side of a region.
- **Merge** (`--merge`, `--merge-distance`): merge overlapping or book-ended
  regions on the same chromosome. `--merge` implies `--sort`.
  `--merge-distance` sets the maximum gap for merging (default: 0).
- **Interval subtraction** (`--subtract`): remove portions of input regions
  that overlap with a BED file, potentially splitting regions.
- **Interval intersection** (`--intersect`): keep only portions of input
  regions that overlap with a BED file.
- **K-mer tiling** (`--tile`, `--step`): tile each region into fixed-width
  windows. `--step` controls stride (default: same as `--tile`).
- **Reverse complement** (`--rc` / `--reverse-complement`): reverse-complement
  all extracted sequences with full IUPAC ambiguity code support. XOR'd with
  BED strand column.
- **Deduplication** (`--dedup`): remove regions with identical chr, start, end
  before extraction.
- **Sorting** (`--sort`): order regions by natural chromosome order then start
  position.
- **Quiet mode** (`-q` / `--quiet`): suppress progress messages, warnings, and
  the progress bar on stderr.
- **Delimiter override** (`--delimiter`): override auto-detection for `--table`
  and `--sv-table`. Accepts `tab`, `comma`, or a single character.
- **FASTQ output** (`--fastq`, `--qual`): emit FASTQ with a constant quality
  character (default: `I`, phred 40).
- **TSV output** (`--tab-out`): emit tab-separated output with columns chr,
  start, end, name, sequence.
- **Statistics mode** (`--stats`): print per-region statistics (TSV) instead of
  extracting sequences. Columns: chr, start, end, name, length, gc_percent,
  n_count, masked_count.
- **Custom headers** (`--name-template`): format FASTA/FASTQ headers with
  placeholders (`{chr}`, `{start}`, `{end}`, `{name}`, `{length}`, `{index}`,
  `{strand}`).
- **FASTA wrapping control** (`--line-width`, `--no-wrap`): configure output
  line width or disable wrapping entirely.
- **Sequence transforms**: `--to-rna` (T to U), `--translate` (standard genetic
  code, stop codons as `*`), `--uppercase`, `--lowercase`.
- **Masking** (`--mask-bed`, `--hard-mask`, `--soft-mask`): mask bases within
  extracted sequences using a BED file. Hard mask (N) is the default;
  `--soft-mask` lowercases instead.
- **No-index mode** (`--no-index`): load the FASTA into memory without an FAI
  index. Required for stdin input (`-`).
- **Multiple FASTA files**: pass several FASTA files as positional arguments.
  Contigs are looked up across all files (each must appear in exactly one).
- **Bgzip/gzip support**: compressed FASTA files are decompressed transparently.
  The FAI index is looked up next to the original file first.
- **Progress bar**: shown on stderr when writing to a file output.
- **Stdin support for `--list`**: pass `-` as the path.
- **Bulk-read optimisation**: nearby regions on the same contig are read in a
  single I/O operation, reducing seek overhead.
- **Streaming ordered output**: stdout and single-file output buffer results in
  memory to preserve input order while extracting in parallel.
- **Library crate** with public re-exports: `Region`, `FaiRecord`,
  `extract_region`, `reverse_complement`, `build_fai`, `read_fai`,
  `parse_region_str`, `wrap_fasta`, `deduplicate_regions`, `sort_regions`,
  `merge_regions`, `tile_regions`, `resolve_flanks`, `parse_regions_bed`,
  `parse_regions_gff`, `parse_regions_vcf`, `intersect_regions`,
  `subtract_regions`, `MaskIndex`, `MaskMode`, `TransformConfig`,
  `expand_template`, `RegionStats`, `compute_stats`, `VcfRecord`, `is_gzip`,
  `resolve_bgzip`.
- **Module split**: source code reorganised into `fai`, `region`, `table`,
  `extract`, `output`, `validate`, `vcf`, `gff`, `intervals`, `mask`,
  `transform`, `template`, `stats`, and `noindex` modules.
- **GitHub Actions CI pipeline** (fmt, clippy, test, security audit) with
  `release/**` branch triggers.
- **MSRV** (1.87) documented in `Cargo.toml` (`rust-version`) and README.
- **CHANGELOG.md** following Keep a Changelog format.
- **BED validation**: rejects empty intervals (start == end) and malformed
  coordinates (start > end).
- **CLI guards**: `--dedup`/`--sort`/`--merge`/`--tile` forbidden with
  `--sv-table --output-dir` to protect paired region ordering. `--delimiter`
  requires `--table` or `--sv-table`. `--flank-left` and `--flank-right` must
  be specified together. `--tile` and `--step` must be > 0.
- **Detection of inconsistent line widths** in `build_fai`.
- **Warning when thread pool construction fails** (falls back to default pool).
- **Clap `conflicts_with` constraints** between mutually exclusive flags.
- **Test expansion**: 31 tests in v0.1.0 grew to 123+ tests (59 unit, 64
  integration) covering BED, VCF, GFF, intervals, tiling, masking, transforms,
  reverse complement, dedup, sort, merge, delimiter sniffing, edge cases, and
  error paths.

### Changed

- `reverse_complement` rewritten to use byte iteration instead of char
  iteration for better performance on ASCII-only sequences.
- File handles reused via thread-local storage across parallel extraction
  instead of opening per region.
- Coordinate semantics clarified in usage documentation (1-based inclusive).
- Code formatting fixed across the codebase.
- Removed `#[allow(dead_code)]` annotations: all public items are now
  exercised from main or re-exported.
- Strand column now recognised in BED (column 6) and table inputs (`STRAND`
  column).

### Fixed

- Division-by-zero guard on malformed FAI entries (`line_bases = 0`).
- `--sv-table` regions now correctly merge with other region sources.
- Flank overflow on upper bound clamped to sequence length via
  `validate_and_clamp_regions`.
- I/O errors in `parse_regions_list` now propagate instead of panicking.
- Removed bogus bgzip support claim from earlier documentation (now properly
  implemented).

## [0.1.0] - 2025-03-25

### Added

- Initial release.
- Multi-sequence extraction from FASTA files using `.fai` indexing.
- Parallel extraction with configurable thread count via Rayon.
- Inline `--regions` flag for comma-separated region specifications.
- `--list` flag for file-based region input (one region per line).
- `--table` flag for CSV/TSV tables with named columns (range and position
  mode).
- `--sv-table` flag for structural variant breakpoint tables.
- `--flank` flag for position-mode tables.
- `--output` and `--output-dir` for single-file and per-region output.
- Automatic `.fai` index building when missing (unless `--no-build-fai`).
- Case-insensitive column name matching.
- Optional `NAME` column for custom output naming.

[0.2.1]: https://github.com/trentzz/multiseqex/compare/v0.2.0...v0.2.1
[0.2.0]: https://github.com/trentzz/multiseqex/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/trentzz/multiseqex/releases/tag/v0.1.0
