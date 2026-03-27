# T001-007: Test FAI build round-trip accuracy

**Epic**: TEST-001
**Priority**: medium
**Depends on**: none
**Status**: done

## Goal

No test verifies that a FAI built by `build_fai` produces correct extraction
results. Create a test that writes a multi-contig FASTA with known content,
auto-builds the FAI, extracts specific ranges, and verifies the extracted
sequences match the expected bases.

## Success Criteria

- [x] Test writes a FASTA with at least two contigs and varying line widths.
- [x] Test extracts ranges spanning line boundaries.
- [x] Extracted sequences are verified against expected content.
- [x] All tests pass.
- [x] /update has been run after changes.
