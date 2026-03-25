# T001-005: Test filename collision in --output-dir

**Epic**: TEST-001
**Priority**: medium
**Depends on**: none
**Status**: done

## Goal

When two regions produce the same filename (e.g. same chr/start/end, or same
NAME/start/end), the second silently overwrites the first. Verify this
behaviour is at least documented by a test, and consider whether a warning or
error is more appropriate.

## Success Criteria

- [x] Integration test demonstrates the collision scenario.
- [x] Behaviour is either: (a) tested and documented, or (b) changed to warn/error.
- [x] All existing tests pass.
- [x] /update has been run after changes.
