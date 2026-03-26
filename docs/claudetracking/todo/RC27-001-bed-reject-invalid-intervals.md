# RC27-001: BED parser should reject start > end and handle start==end

**Epic**: REVIEW-CYCLE-2026-03-27
**Priority**: high
**Depends on**: none
**Status**: todo

## Goal

The BED parser currently accepts lines where `start > end`, which produces
nonsensical regions. Add validation to reject these with a clear error
message. Also handle the `start == end` case (an empty interval in BED
coordinates). Either skip empty intervals with a warning or reject them,
depending on what makes sense for extraction.

## Success Criteria

- [ ] BED lines with `start > end` produce a clear error and stop parsing.
- [ ] BED lines with `start == end` are handled explicitly (skipped with
      warning or rejected with error).
- [ ] Existing valid BED files still parse correctly.
- [ ] Unit tests cover both the `start > end` and `start == end` cases.
- [ ] All tests pass.
- [ ] `/update` has been run after changes.
