# F001-003: Add progress reporting for large extractions

**Epic**: FEAT-001
**Priority**: low
**Depends on**: none
**Status**: todo

## Goal

For large region sets, the tool gives no feedback during extraction. Add a
simple progress indicator to stderr showing the number of regions processed
out of the total. Use a lightweight approach (e.g. atomic counter printed
every N regions) to avoid performance overhead.

## Success Criteria

- [ ] Progress messages appear on stderr for extractions with more than 100 regions.
- [ ] Progress does not appear for small extractions (avoids noise).
- [ ] A `--quiet` flag suppresses progress output.
- [ ] No measurable performance regression.
- [ ] All existing tests pass.
- [ ] /update has been run after changes.
