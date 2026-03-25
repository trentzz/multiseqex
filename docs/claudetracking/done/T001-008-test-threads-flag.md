# T001-008: Test --threads flag

**Epic**: TEST-001
**Priority**: low
**Depends on**: none
**Status**: done

## Goal

The `--threads` flag is untested. Add an integration test that passes
`--threads 1` and verifies correct output. This confirms single-threaded mode
works and that the flag is accepted.

## Success Criteria

- [x] Integration test passes `--threads 1` with a multi-region extraction.
- [x] Output is correct.
- [x] All tests pass.
- [x] /update has been run after changes.
