# RC25-012: Clarify coordinate semantics in docs

**Epic**: REVIEW-CYCLE-2026-03-25
**Priority**: medium
**Depends on**: none
**Status**: done

## Goal

Documentation does not state the coordinate system (1-based inclusive). The
flank example in `docs/usage.md` is correct but ambiguous without stating the
convention. Add explicit coordinate semantics to the usage docs and README.

## Success Criteria

- [x] Coordinate system (1-based inclusive) is stated in usage docs.
- [x] The flank example includes a note about the semantics.
- [x] All tests pass.
- [x] `/update` has been run after changes.
