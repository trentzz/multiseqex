# F002-002: Deduplicate identical regions

**Epic**: FEAT-002
**Priority**: low
**Depends on**: none
**Status**: done

## Goal

When the same region appears multiple times (e.g. from combining `--regions`
and `--table`, or from a table with duplicate rows), the tool extracts and
outputs it multiple times. Add a `--dedup` flag that removes exact duplicates
and warns about them on stderr.

## Success Criteria

- [x] `--dedup` flag removes regions with identical chr, start, end.
- [x] Duplicate count is reported on stderr.
- [x] Without `--dedup`, behaviour is unchanged (duplicates preserved).
- [x] Integration test covers deduplication.
- [x] All tests pass.
- [x] /update has been run after changes.
