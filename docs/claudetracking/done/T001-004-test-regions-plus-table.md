# T001-004: Test combining --regions with --table

**Epic**: TEST-001
**Priority**: medium
**Depends on**: none
**Status**: done

## Goal

An integration test exists for `--regions` + `--list`, but not for `--regions`
+ `--table`. Add one to verify both sources are merged.

## Success Criteria

- [x] Integration test combines --regions and --table and asserts regions from both appear in output.
- [x] All existing tests pass.
- [x] /update has been run after changes.
