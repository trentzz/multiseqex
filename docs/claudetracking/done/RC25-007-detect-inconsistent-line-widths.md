# RC25-007: Detect inconsistent line widths in build_fai

**Epic**: REVIEW-CYCLE-2026-03-25
**Priority**: medium
**Depends on**: none
**Status**: done

## Goal

`build_fai` records `line_bases`/`line_bytes` from only the first sequence line.
If a FASTA has inconsistent line widths, the FAI is wrong and extraction silently
returns incorrect sequences. Detect inconsistency and warn or error.

## Success Criteria

- [x] `build_fai` checks that all non-final sequence lines have the same width.
- [x] A warning or error is emitted for inconsistent widths.
- [x] A test covers the inconsistent-width case.
- [x] All tests pass.
- [x] `/update` has been run after changes.
