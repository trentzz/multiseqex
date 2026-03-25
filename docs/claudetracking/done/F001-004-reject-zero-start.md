# F001-004: Reject zero-start coordinates in 1-based mode

**Epic**: FEAT-001
**Priority**: high
**Depends on**: none
**Status**: done

## Goal

The coordinate system is documented as 1-based inclusive. Passing `chr1:0-10`
is silently accepted. In `extract_region`, `pos` starts at `r.start` and the
byte offset calculation uses `(pos - 1) / lb`. When `pos = 0`, this underflows
(wraps to `u64::MAX`), producing a seek to an invalid offset.

Add validation in `parse_region_str` and table parsing to reject start or end
values of 0 with a clear error message.

## Success Criteria

- [x] `parse_region_str("chr1:0-10", None)` returns an error mentioning 1-based coordinates.
- [x] Table rows with START=0 or END=0 produce an error.
- [x] Unit tests cover both cases.
- [x] All existing tests pass.
- [x] /update has been run after changes.
