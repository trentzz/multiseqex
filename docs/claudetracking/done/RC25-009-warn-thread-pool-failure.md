# RC25-009: Warn on thread pool build failure

**Epic**: REVIEW-CYCLE-2026-03-25
**Priority**: medium
**Depends on**: none
**Status**: done

## Goal

Thread pool initialisation failure is silently discarded with `.ok()`. If the
user requests `--threads 4` and it fails, they get the default with no warning.
Print a warning on failure.

## Success Criteria

- [x] A warning is printed to stderr if the thread pool fails to build.
- [x] All tests pass.
- [x] `/update` has been run after changes.
