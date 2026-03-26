# RC27-012: CI triggers and code hygiene batch

**Epic**: REVIEW-CYCLE-2026-03-27
**Priority**: low
**Depends on**: none
**Status**: done

## Goal

Address three low-priority housekeeping items:

- L1: Remove `#[allow(dead_code)]` annotations by adding doc-tests that
  exercise the annotated items, proving they are part of the public API.
- L2: Add a comment to the natural sort implementation explaining the
  overflow risk on extremely large numeric substrings, and why it is
  acceptable in practice.
- L3: Add `release/**` to the CI push triggers so that release branches
  get the same checks as `main` and `dev`.

## Success Criteria

- [x] No `#[allow(dead_code)]` annotations remain (items either have
      doc-tests or are genuinely removed).
- [x] Natural sort has a comment explaining the overflow edge case.
- [x] CI workflow triggers on `release/**` branches.
- [x] All tests pass.
- [x] `/update` has been run after changes.
