# ROBUST-001: Robustness and Edge-Case Hardening

## Goal

Harden the tool against malformed input, edge cases, and silent data corruption.

## Motivation

Several edge cases can produce silently wrong results or confusing errors:
duplicate contig names in FASTA, empty header lines, stale FAI indexes, and
inconsistent test fixtures. Fixing these improves reliability for production
bioinformatics pipelines.

## Scope

- Fix test fixture with inconsistent line widths.
- Handle duplicate contig names in FASTA/FAI.
- Handle empty contig names (bare `>` header).
- Detect stale FAI indexes (mtime check).
- Guard `extract_region` against zero-start input.

## Tasks

- R001-001: Fix chr2 inconsistent line widths in test fixture
- R001-002: Warn on duplicate contig names in FAI
- R001-003: Handle empty contig names in FASTA header
- R001-004: Detect stale FAI index
- R001-005: Guard extract_region against zero-start
