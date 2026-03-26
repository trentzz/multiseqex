# RC27-011: Remove dead Strand abstraction or wire it up

**Epic**: REVIEW-CYCLE-2026-03-27
**Priority**: low
**Depends on**: none
**Status**: todo

## Goal

The codebase contains a `Strand` type that is defined but never used in the
extraction pipeline. Either remove it entirely to reduce dead code, or wire
it up properly: parse the BED strand column (column 6), store strand on
`Region`, and apply per-region reverse complement for minus-strand regions.

Removing is simpler. Wiring up is more useful but larger in scope.

## Success Criteria

- [ ] `Strand` is either removed or fully integrated into BED parsing and
      extraction.
- [ ] No `#[allow(dead_code)]` annotations remain for strand-related code.
- [ ] If wired up: minus-strand BED regions are automatically reverse
      complemented.
- [ ] If removed: no references to `Strand` remain in the codebase.
- [ ] All tests pass.
- [ ] `/update` has been run after changes.
