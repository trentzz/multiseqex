# A001-002: Extract library crate (dual crate)

**Epic**: ARCH-001
**Priority**: medium
**Depends on**: A001-001
**Status**: todo

## Goal

Restructure the project as a dual crate: `src/lib.rs` exposes the core API
(FAI parsing, region parsing, sequence extraction) and `src/main.rs` is a thin
CLI wrapper. This enables programmatic use and cleaner testing.

## Success Criteria

- [ ] `src/lib.rs` exists and exports core types and functions.
- [ ] `src/main.rs` imports from the library and contains only CLI logic.
- [ ] All existing tests pass.
- [ ] `cargo doc` generates documentation for the library API.
- [ ] `/update` has been run after changes.
