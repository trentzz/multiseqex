# REVIEW-CYCLE-2026-03-25

## Goal

Critical review of multiseqex v0.1.0. Fix correctness bugs, improve robustness,
close test gaps, and clean up documentation.

## Motivation

First formal review cycle after initial release. The Opus review found 5 high-priority
issues (including a false bgzip claim, division-by-zero risk, and file descriptor
exhaustion), 8 medium-priority issues (zero unit tests, silent error swallowing,
overflow bugs), and 8 low-priority polish items.

## Scope

All findings from the 2026-03-25 review. High and medium priority items are
in-scope. Low priority items are included where they can be fixed alongside
related work.

## Tasks

- RC25-001: Remove bogus bgzip claim (H1)
- RC25-002: Guard against division-by-zero on malformed FAI (H2)
- RC25-003: Reuse file handles in parallel extraction (H3)
- RC25-004: Fix --sv-table mixing with other region sources (H4)
- RC25-005: Fix cargo fmt failures (H5)
- RC25-006: Add unit tests for core functions (M1)
- RC25-007: Detect inconsistent line widths in build_fai (M2)
- RC25-008: Propagate I/O errors in parse_regions_list (M3)
- RC25-009: Warn on thread pool build failure (M4)
- RC25-010: Use clap conflicts_with for output args (M5)
- RC25-011: Fix flank overflow on upper bound (M6)
- RC25-012: Clarify coordinate semantics in docs (M7)
- RC25-013: Low-priority polish batch (L1-L8)
