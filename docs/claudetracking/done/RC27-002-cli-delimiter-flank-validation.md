# RC27-002: Add CLI validation for --delimiter and --flank

**Epic**: REVIEW-CYCLE-2026-03-27
**Priority**: high
**Depends on**: none
**Status**: done

## Goal

Add a `requires` constraint so that `--delimiter` is only accepted when
`--table` or `--sv-table` is also provided. Without a table mode, the
delimiter has no effect and silently misleads.

Also validate `--flank` usage. If `--flank` is provided without position-mode
regions (i.e. only name-mode regions from `--list` or `--table`), warn or
error. Flanking only applies to regions with numeric coordinates.

## Success Criteria

- [x] `--delimiter` without `--table` or `--sv-table` produces a clear error.
- [x] `--flank` without any position-mode regions produces a warning or error.
- [x] Valid flag combinations still work correctly.
- [x] Tests cover the rejected combinations.
- [x] All tests pass.
- [x] `/update` has been run after changes.
