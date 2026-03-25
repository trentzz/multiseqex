# RC25-005: Fix cargo fmt failures

**Epic**: REVIEW-CYCLE-2026-03-25
**Priority**: high
**Depends on**: none
**Status**: done

## Goal

`cargo fmt -- --check` fails. Run `cargo fmt` to fix formatting in `src/main.rs`
and `tests/integration.rs`.

## Success Criteria

- [x] `cargo fmt -- --check` passes.
- [x] All tests pass.
- [x] `/update` has been run after changes.
