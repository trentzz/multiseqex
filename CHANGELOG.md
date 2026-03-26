# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0] - 2026-03-27

### Added

- BED file input (`--bed`) with 0-based half-open to 1-based inclusive
  coordinate conversion. Supports optional name column, comment lines,
  and blank lines.
- BED flanking support: `--flank` extends BED regions symmetrically after
  coordinate conversion.
- Reverse complement flag (`--rc` / `--reverse-complement`) with full IUPAC
  ambiguity code support.
- Deduplication flag (`--dedup`) removes regions with identical chr, start, end
  before extraction.
- Sort flag (`--sort`) orders regions by natural chromosome order then start
  position.
- Quiet mode (`-q` / `--quiet`) suppresses progress messages and warnings on
  stderr.
- Delimiter override (`--delimiter`) for `--table` and `--sv-table`. Accepts
  `tab`, `comma`, or a single character. Auto-detection sniffs for tabs in the
  header line when the file extension is not `.tsv`.
- Stdin support for `--list` (pass `-` as the path).
- Bulk-read optimisation: nearby regions on the same contig are read in a
  single I/O operation, reducing seek overhead.
- Streaming ordered output: stdout and single-file output buffer results in
  memory to preserve input order while extracting in parallel.
- Library crate with public re-exports (`Region`, `FaiRecord`, `extract_region`,
  `reverse_complement`, `build_fai`, `read_fai`, `parse_region_str`,
  `wrap_fasta`, `deduplicate_regions`, `sort_regions`, `parse_regions_bed`).
- Module split: source code reorganised into `fai`, `region`, `table`,
  `extract`, `output`, and `validate` modules.
- GitHub Actions CI pipeline (fmt, clippy, test, security audit) with
  `release/**` branch triggers.
- MSRV (1.87) documented in `Cargo.toml` (`rust-version`) and README.
- CHANGELOG.md following Keep a Changelog format.
- BED validation: rejects empty intervals (start == end) and malformed
  coordinates (start > end).
- CLI guards: `--dedup`/`--sort` forbidden with `--sv-table --output-dir` to
  protect paired region ordering. `--delimiter` requires `--table` or
  `--sv-table`.
- Detection of inconsistent line widths in `build_fai`.
- Warning when thread pool construction fails (falls back to default pool).
- Clap `conflicts_with` constraint between `--output` and `--output-dir`.
- Test expansion: 31 tests in v0.1.0 grew to 123+ tests (59 unit, 64
  integration) covering BED, reverse complement, dedup, sort, delimiter
  sniffing, edge cases, and error paths.

### Changed

- `reverse_complement` rewritten to use byte iteration instead of char
  iteration for better performance on ASCII-only sequences.
- File handles reused via thread-local storage across parallel extraction
  instead of opening per region.
- Coordinate semantics clarified in usage documentation (1-based inclusive).
- Code formatting fixed across the codebase.
- Removed `#[allow(dead_code)]` annotations: all public items are now
  exercised from main or re-exported.
- Strand column removed from region parsing (was unused and misleading).

### Fixed

- Division-by-zero guard on malformed FAI entries (`line_bases = 0`).
- `--sv-table` regions now correctly merge with other region sources.
- Flank overflow on upper bound clamped to sequence length via
  `validate_and_clamp_regions`.
- I/O errors in `parse_regions_list` now propagate instead of panicking.
- Removed bogus bgzip support claim from documentation.

## [0.1.0] - 2025-03-25

### Added

- Initial release.
- Multi-sequence extraction from FASTA files using `.fai` indexing.
- Parallel extraction with configurable thread count via Rayon.
- Inline `--regions` flag for comma-separated region specifications.
- `--list` flag for file-based region input (one region per line).
- `--table` flag for CSV/TSV tables with named columns (range and position mode).
- `--sv-table` flag for structural variant breakpoint tables.
- `--flank` flag for position-mode tables.
- `--output` and `--output-dir` for single-file and per-region output.
- Automatic `.fai` index building when missing (unless `--no-build-fai`).
- Case-insensitive column name matching.
- Optional `NAME` column for custom output naming.

[0.2.0]: https://github.com/trentzz/multiseqex/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/trentzz/multiseqex/releases/tag/v0.1.0
