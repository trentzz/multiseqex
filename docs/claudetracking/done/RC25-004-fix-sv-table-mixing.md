# RC25-004: Fix --sv-table mixing with other region sources

**Epic**: REVIEW-CYCLE-2026-03-25
**Priority**: high
**Depends on**: none
**Status**: done

## Goal

When `--sv-table` is combined with `--regions` or `--table`, the `is_sv` flag
is set to true for all regions. Non-SV regions get mixed into the SV vector,
corrupting `--output-dir` output (which expects paired regions). Either make
`--sv-table` conflict with other region sources, or track SV vs non-SV regions
separately.

## Success Criteria

- [x] `--sv-table` conflicts with `--regions` and `--table`, OR SV/non-SV
      regions are tracked and output separately.
- [x] A test covers the conflict or the separate handling.
- [x] All tests pass.
- [x] `/update` has been run after changes.
