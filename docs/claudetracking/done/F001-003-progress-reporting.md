# F001-003: Add progress reporting for large extractions

**Epic**: FEAT-001
**Priority**: low
**Depends on**: none
**Status**: done

## Goal

For large region sets, the tool gives no feedback during extraction. Add a
simple progress indicator to stderr showing the number of regions processed
out of the total. Use a lightweight approach (e.g. atomic counter printed
every N regions) to avoid performance overhead.

## Success Criteria

- [x] Progress messages appear on stderr for extractions with more than 100 regions.
- [x] Progress does not appear for small extractions (avoids noise).
- [x] A `--quiet` flag suppresses progress output.
- [x] No measurable performance regression.
- [x] All existing tests pass.
- [x] /update has been run after changes.
