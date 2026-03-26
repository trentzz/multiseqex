# RC27-002: Add CLI validation for --delimiter and --flank

**Epic**: REVIEW-CYCLE-2026-03-27
**Priority**: high
**Depends on**: none
**Status**: todo

## Goal

Add a `requires` constraint so that `--delimiter` is only accepted when
`--table` or `--sv-table` is also provided. Without a table mode, the
delimiter has no effect and silently misleads.

Also validate `--flank` usage. If `--flank` is provided without position-mode
regions (i.e. only name-mode regions from `--list` or `--table`), warn or
error. Flanking only applies to regions with numeric coordinates.

## Success Criteria

- [ ] `--delimiter` without `--table` or `--sv-table` produces a clear error.
- [ ] `--flank` without any position-mode regions produces a warning or error.
- [ ] Valid flag combinations still work correctly.
- [ ] Tests cover the rejected combinations.
- [ ] All tests pass.
- [ ] `/update` has been run after changes.
