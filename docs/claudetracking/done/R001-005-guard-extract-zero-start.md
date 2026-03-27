# R001-005: Guard extract_region against zero-start

**Epic**: ROBUST-001
**Priority**: high
**Depends on**: none
**Status**: done

## Goal

`extract_region` computes `(pos - 1)` on line 344 without guarding against
`r.start == 0`. While `validate_and_clamp_regions` clamps start to 1, the
extraction function itself should be defensively safe. A debug assertion or
early error protects against future callers that bypass validation.

## Success Criteria

- [x] `extract_region` returns an error (not a panic) if `r.start == 0`.
- [x] A unit test covers the zero-start guard.
- [x] All tests pass.
- [x] /update has been run after changes.
