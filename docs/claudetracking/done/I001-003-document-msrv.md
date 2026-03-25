# I001-003: Document MSRV (Rust 1.87+)

**Epic**: INFRA-001
**Priority**: low
**Depends on**: none
**Status**: done

## Goal

The project uses `edition = "2024"` which requires Rust 1.87+. This is not
documented in Cargo.toml (`rust-version` field) or prominently in the README.

## Success Criteria

- [x] `rust-version = "1.87"` added to Cargo.toml `[package]`.
- [x] README prerequisites section mentions Rust 1.87+.
- [x] /update has been run after changes.
