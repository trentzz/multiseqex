# RC27-010: Add BED integration tests

**Epic**: REVIEW-CYCLE-2026-03-27
**Priority**: low
**Depends on**: RC27-001
**Status**: done

## Goal

The BED input path lacks integration tests that exercise the full pipeline
from BED file to extracted sequences. Add tests covering: a well-formed BED
file, a BED file with extra columns (ignored), edge cases (first base, last
base of a contig), and interaction with other flags (`--rc`, `--flank`).

## Success Criteria

- [x] At least 4 integration tests exercise `--bed` end-to-end.
- [x] Tests cover: basic extraction, extra columns, edge coordinates, and
      flag interaction.
- [x] All tests pass.
- [x] `/update` has been run after changes.
