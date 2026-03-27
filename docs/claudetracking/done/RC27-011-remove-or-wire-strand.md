# RC27-011: Remove dead Strand abstraction or wire it up

**Epic**: REVIEW-CYCLE-2026-03-27
**Priority**: low
**Depends on**: none
**Status**: done

## Goal

The codebase contains a `Strand` type that is defined but never used in the
extraction pipeline. Either remove it entirely to reduce dead code, or wire
it up properly: parse the BED strand column (column 6), store strand on
`Region`, and apply per-region reverse complement for minus-strand regions.

Removing is simpler. Wiring up is more useful but larger in scope.

## Success Criteria

- [x] `Strand` is either removed or fully integrated into BED parsing and
      extraction.
- [x] No `#[allow(dead_code)]` annotations remain for strand-related code.
- [x] If wired up: minus-strand BED regions are automatically reverse
      complemented.
- [x] If removed: no references to `Strand` remain in the codebase.
- [x] All tests pass.
- [x] `/update` has been run after changes.
