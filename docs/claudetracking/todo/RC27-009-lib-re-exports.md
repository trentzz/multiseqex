# RC27-009: Add top-level re-exports in lib.rs

**Epic**: REVIEW-CYCLE-2026-03-27
**Priority**: medium
**Depends on**: none
**Status**: todo

## Goal

Users of the library crate must currently navigate module paths to reach
common types. Add top-level `pub use` re-exports in `lib.rs` for the most
commonly used types and functions so that `use multiseqex::Region` works
without knowing the module layout.

## Success Criteria

- [ ] Common types (e.g. `Region`, `Fai`, key parse/extract functions) are
      re-exported from the crate root.
- [ ] Existing module paths still work (re-exports are additive).
- [ ] Doc-tests or examples use the short import paths.
- [ ] All tests pass.
- [ ] `/update` has been run after changes.
