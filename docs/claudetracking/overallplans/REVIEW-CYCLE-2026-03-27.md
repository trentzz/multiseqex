# REVIEW-CYCLE-2026-03-27

## Goal

Second review cycle for multiseqex. Fix validation gaps, improve CLI
correctness, address streaming limitations, and update documentation for
v0.2.0 features.

## Motivation

Opus review of the codebase after the v0.1.0 release cycle found 4
high-priority issues (BED parser accepts invalid intervals, missing CLI
argument validation, incorrect dedup/sort behaviour with SV pairing, and
memory buffering disguised as streaming output). Medium-priority findings
cover documentation staleness, missing BED flanking support, a suboptimal
reverse complement path, and API ergonomics. Low-priority items address dead
code and CI configuration.

## Scope

All findings from the 2026-03-27 review. High and medium priority items are
in-scope. Low priority items are included where they can be fixed alongside
related work.

## Tasks

- RC27-001: BED parser should reject start > end and handle start==end (H1, M4)
- RC27-002: Add CLI validation for --delimiter and --flank (H4, M1)
- RC27-003: Forbid --dedup/--sort with --sv-table --output-dir (M2)
- RC27-004: Fix streaming output to be truly streaming (H3)
- RC27-005: Update README.md and usage docs for v0.2.0 flags (M6)
- RC27-006: Update CHANGELOG.md with v0.2.0 features (M7)
- RC27-007: Add BED flanking support (M3)
- RC27-008: Optimise reverse_complement byte-level iteration (M5)
- RC27-009: Add top-level re-exports in lib.rs (M8)
- RC27-010: Add BED integration tests (L8)
- RC27-011: Remove dead Strand abstraction or wire it up (L6)
- RC27-012: CI triggers and code hygiene batch (L1, L2, L3)
