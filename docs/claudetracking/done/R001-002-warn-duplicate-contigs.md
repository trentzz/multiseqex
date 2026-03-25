# R001-002: Warn on duplicate contig names in FAI

**Epic**: ROBUST-001
**Priority**: medium
**Depends on**: none
**Status**: done

## Goal

`read_fai` uses a HashMap keyed by contig name. If a FASTA has duplicate contig
names, the last entry silently overwrites earlier ones. This produces wrong
results with no warning. Detect duplicates in `read_fai` and in `build_fai`,
and emit a warning or error.

## Success Criteria

- [x] `read_fai` warns or errors when duplicate contig names are encountered.
- [x] `build_fai` warns or errors when duplicate contig names are encountered.
- [x] A unit test covers the duplicate contig case.
- [x] All tests pass.
- [x] /update has been run after changes.
