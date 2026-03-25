# R001-003: Handle empty contig names in FASTA header

**Epic**: ROBUST-001
**Priority**: low
**Depends on**: none
**Status**: done

## Goal

A FASTA header line containing only `>` (no sequence name) causes
`parse_fasta_header` to return an empty string. This creates a FAI entry with
an empty key, which is technically valid but likely indicates a malformed file.
Emit a warning during FAI building when an empty contig name is encountered.

## Success Criteria

- [x] `build_fai` warns when a contig has an empty name.
- [x] A unit test covers the bare `>` header case.
- [x] All tests pass.
- [x] /update has been run after changes.
