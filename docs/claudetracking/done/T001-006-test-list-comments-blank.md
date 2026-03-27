# T001-006: Test --list with comments and blank lines

**Epic**: TEST-001
**Priority**: medium
**Depends on**: none
**Status**: done

## Goal

The `parse_regions_list` function skips blank lines and lines starting with `#`,
but no integration test verifies this behaviour. Add a test using a list file
with comments, blank lines, and valid regions interleaved.

## Success Criteria

- [x] Integration test uses a list file with `#` comments and blank lines.
- [x] Test verifies only valid regions are extracted.
- [x] All tests pass.
- [x] /update has been run after changes.
