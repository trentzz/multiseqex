# RC25-013: Low-priority polish batch

**Epic**: REVIEW-CYCLE-2026-03-25
**Priority**: low
**Depends on**: none
**Status**: done

## Goal

Address low-priority findings from the review:

- L1: Fix `documentation` URL in Cargo.toml (point to docs.rs or remove).
- L4: Guard `wrap_fasta` against `width == 0`.
- L5: Fix misleading comment on `count_bases` (says IUPAC but matches all
  ASCII alphabetic).
- L8: Align author email in Cargo.toml with git config.

L2 (deny missing_docs), L3 (Region clone), L6 (version test), L7 (delimiter
sniffing) are deferred as too low impact.

## Success Criteria

- [x] `documentation` field in Cargo.toml is correct or removed.
- [x] `wrap_fasta` handles `width == 0` without panic.
- [x] `count_bases` comment accurately describes behaviour.
- [x] Author email is consistent.
- [x] All tests pass.
- [x] `/update` has been run after changes.
