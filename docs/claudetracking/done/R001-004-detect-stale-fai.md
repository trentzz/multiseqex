# R001-004: Detect stale FAI index

**Epic**: ROBUST-001
**Priority**: medium
**Depends on**: none
**Status**: done

## Goal

If the FASTA file is modified after its FAI was built, the index is stale and
extractions produce silently wrong results. Compare mtime of the FASTA and FAI
files. If the FASTA is newer, warn the user and optionally rebuild.

## Success Criteria

- [x] When the FASTA file mtime is newer than the FAI mtime, a warning is
      printed to stderr.
- [x] If `--no-build-fai` is not set, the stale FAI is rebuilt automatically
      (with a message).
- [x] A test covers the stale-index detection path.
- [x] All tests pass.
- [x] /update has been run after changes.
