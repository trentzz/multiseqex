# F001-001: Add BED file input support

**Epic**: FEAT-001
**Priority**: medium
**Depends on**: none
**Status**: todo

## Goal

BED is the most common genomic interval format. It uses 0-based, half-open
coordinates (start is 0-based inclusive, end is exclusive). Add a `--bed`
flag that reads a BED file and converts to the internal 1-based inclusive
representation.

Support BED3 (chrom, start, end) at minimum. If a fourth column (name) is
present, use it as the region name.

## Success Criteria

- [ ] `--bed` flag accepts a BED file path.
- [ ] BED coordinates are correctly converted to 1-based inclusive (start+1, end unchanged).
- [ ] BED4 name column is used as region name when present.
- [ ] Integration tests cover BED3 and BED4 input.
- [ ] Documentation updated (README and usage.md).
- [ ] All existing tests pass.
- [ ] /update has been run after changes.
