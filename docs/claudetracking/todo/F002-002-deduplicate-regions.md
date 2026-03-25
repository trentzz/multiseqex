# F002-002: Deduplicate identical regions

**Epic**: FEAT-002
**Priority**: low
**Depends on**: none
**Status**: todo

## Goal

When the same region appears multiple times (e.g. from combining `--regions`
and `--table`, or from a table with duplicate rows), the tool extracts and
outputs it multiple times. Add a `--dedup` flag that removes exact duplicates
and warns about them on stderr.

## Success Criteria

- [ ] `--dedup` flag removes regions with identical chr, start, end.
- [ ] Duplicate count is reported on stderr.
- [ ] Without `--dedup`, behaviour is unchanged (duplicates preserved).
- [ ] Integration test covers deduplication.
- [ ] All tests pass.
- [ ] /update has been run after changes.
