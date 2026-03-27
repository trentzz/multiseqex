# T001-002: Test malformed CSV rows

**Epic**: TEST-001
**Priority**: medium
**Depends on**: none
**Status**: done

## Goal

Add integration tests for malformed table input: rows with missing fields,
non-numeric values in numeric columns, and rows with extra fields. Verify the
tool produces a clear error rather than panicking.

## Success Criteria

- [x] Test: table row with missing END field produces an error.
- [x] Test: table row with non-numeric START produces an error mentioning the row.
- [x] Test: extra fields in a row do not cause errors.
- [x] All existing tests pass.
- [x] /update has been run after changes.
