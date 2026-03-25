# RC25-011: Fix flank overflow on upper bound

**Epic**: REVIEW-CYCLE-2026-03-25
**Priority**: medium
**Depends on**: none
**Status**: done

## Goal

`end: pos + flank` uses plain addition which can overflow `u64`. The lower bound
uses `saturating_sub` but the upper bound does not use `saturating_add`. Fix all
occurrences to use `saturating_add`.

## Success Criteria

- [x] All `pos + flank` and similar upper-bound calculations use `saturating_add`.
- [x] A test covers a very large position + flank value.
- [x] All tests pass.
- [x] `/update` has been run after changes.
