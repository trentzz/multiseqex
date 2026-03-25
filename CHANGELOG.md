# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0] - 2026-03-26

### Added

- Unit tests for core functions (FAI parsing, region parsing, extraction).
- Detection of inconsistent line widths in `build_fai`.
- Warning when thread pool construction fails (falls back to default pool).
- Clap `conflicts_with` constraint between `--output` and `--output-dir`.
- GitHub Actions CI pipeline (fmt, clippy, test, audit).
- CHANGELOG.md following Keep a Changelog format.
- MSRV (1.87) documented in Cargo.toml (`rust-version`) and README.

### Fixed

- Division-by-zero guard on malformed FAI entries.
- `--sv-table` regions now correctly merge with other region sources.
- Flank overflow on upper bound clamped to sequence length.
- I/O errors in `parse_regions_list` now propagate instead of panicking.
- Removed bogus bgzip support claim from documentation.

### Changed

- File handles reused across parallel extraction instead of opening per region.
- Coordinate semantics clarified in usage documentation (1-based inclusive).
- Code formatting fixed across the codebase.
- Low-priority polish: tidied warnings, naming, and documentation.

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
