# F002-003: Add sorted output mode

**Epic**: FEAT-002
**Priority**: low
**Depends on**: none
**Status**: done

## Goal

Add a `--sort` flag that sorts output regions by genomic coordinate (chr then
start position). This simplifies downstream processing and makes output
deterministic regardless of input order. Use natural sort order for chromosome
names (chr1, chr2, ..., chr10 rather than chr1, chr10, chr2).

## Success Criteria

- [x] `--sort` flag sorts output by chromosome (natural order) then start.
- [x] Without `--sort`, output order matches input order.
- [x] Integration test verifies sorted output.
- [x] All tests pass.
- [x] /update has been run after changes.
