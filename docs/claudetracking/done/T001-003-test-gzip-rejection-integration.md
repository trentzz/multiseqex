# T001-003: Integration test for gzip rejection

**Epic**: TEST-001
**Priority**: medium
**Depends on**: none
**Status**: done

## Goal

The gzip detection is unit-tested but not integration-tested. Add a CLI-level
test that passes a gzip-magic file as the FASTA argument and verifies the error
message.

## Success Criteria

- [x] Integration test creates a temp file with gzip magic bytes and asserts failure with "gzip" in stderr.
- [x] All existing tests pass.
- [x] /update has been run after changes.
