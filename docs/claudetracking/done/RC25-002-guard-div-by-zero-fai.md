# RC25-002: Guard against division-by-zero on malformed FAI

**Epic**: REVIEW-CYCLE-2026-03-25
**Priority**: high
**Depends on**: none
**Status**: done

## Goal

`extract_region` divides by `line_bases` (line 306) without checking for zero.
A malformed FAI with `line_bases=0` causes a panic. Add a guard that returns a
clear error. Also clean up the dead guard code at lines 287-295 that can never
trigger after `validate_and_clamp_regions`.

## Success Criteria

- [x] `extract_region` returns an error (not a panic) when `line_bases == 0`.
- [x] Dead guard code after clamping is removed or documented.
- [x] A test covers the `line_bases=0` case.
- [x] All tests pass.
- [x] `/update` has been run after changes.
