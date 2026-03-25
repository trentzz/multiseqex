# F002-003: Add sorted output mode

**Epic**: FEAT-002
**Priority**: low
**Depends on**: none
**Status**: todo

## Goal

Add a `--sort` flag that sorts output regions by genomic coordinate (chr then
start position). This simplifies downstream processing and makes output
deterministic regardless of input order. Use natural sort order for chromosome
names (chr1, chr2, ..., chr10 rather than chr1, chr10, chr2).

## Success Criteria

- [ ] `--sort` flag sorts output by chromosome (natural order) then start.
- [ ] Without `--sort`, output order matches input order.
- [ ] Integration test verifies sorted output.
- [ ] All tests pass.
- [ ] /update has been run after changes.
