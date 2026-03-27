# R001-001: Fix chr2 inconsistent line widths in test fixture

**Epic**: ROBUST-001
**Priority**: medium
**Depends on**: none
**Status**: done

## Goal

The test fixture `tests/fixtures/test.fa` has inconsistent line widths for chr2
(line 1 has 62 bases, line 2 has 63 bases). The pre-built FAI records
line_bases=62, which is only correct for the first line. Extraction happens to
produce correct results for the current test ranges, but this is fragile. Fix
the fixture so all non-final lines have the same width, and regenerate the FAI.

## Success Criteria

- [x] All non-final sequence lines in `tests/fixtures/test.fa` have consistent
      widths within each contig.
- [x] `tests/fixtures/test.fa.fai` matches the corrected FASTA.
- [x] Total bases per contig remain unchanged (tests expect specific lengths).
- [x] All tests pass.
- [x] /update has been run after changes.
