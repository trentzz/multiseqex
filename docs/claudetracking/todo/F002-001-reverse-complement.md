# F002-001: Add reverse complement support

**Epic**: FEAT-002
**Priority**: medium
**Depends on**: A001-001
**Status**: todo

## Goal

Add the ability to extract reverse complement sequences. This is essential for
strand-aware analysis (e.g. extracting promoter regions on the minus strand).
Support via a `--rc` flag for all regions, or a `STRAND` column in table input
(values: `+`, `-`, `.`).

## Success Criteria

- [ ] A `STRAND` column in `--table` and `--sv-table` controls strand.
- [ ] `--rc` flag reverse-complements all extracted sequences.
- [ ] FASTA headers indicate strand (e.g. `>chr1:100-200(-)` for minus).
- [ ] Unit tests cover complement logic for all IUPAC bases.
- [ ] Integration tests cover `--rc` and `STRAND` column.
- [ ] All tests pass.
- [ ] /update has been run after changes.
