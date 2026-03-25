# T001-001: Test empty table (headers only)

**Epic**: TEST-001
**Priority**: medium
**Depends on**: none
**Status**: done

## Goal

Verify behaviour when a CSV/TSV table has a valid header row but no data rows.
Currently this would produce zero regions, which triggers the "No regions
provided" error. This path should be tested explicitly.

## Success Criteria

- [x] Integration test passes an empty-body table via --table and asserts the error message.
- [x] Same for --sv-table.
- [x] All existing tests pass.
- [x] /update has been run after changes.
