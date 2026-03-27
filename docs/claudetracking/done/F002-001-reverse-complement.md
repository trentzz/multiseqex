# F002-001: Add reverse complement support

**Epic**: FEAT-002
**Priority**: medium
**Depends on**: A001-001
**Status**: done

## Goal

Add the ability to extract reverse complement sequences. This is essential for
strand-aware analysis (e.g. extracting promoter regions on the minus strand).
Support via a `--rc` flag for all regions, or a `STRAND` column in table input
(values: `+`, `-`, `.`).

## Success Criteria

- [x] A `STRAND` column in `--table` and `--sv-table` controls strand.
- [x] `--rc` flag reverse-complements all extracted sequences.
- [x] FASTA headers indicate strand (e.g. `>chr1:100-200(-)` for minus).
- [x] Unit tests cover complement logic for all IUPAC bases.
- [x] Integration tests cover `--rc` and `STRAND` column.
- [x] All tests pass.
- [x] /update has been run after changes.
