# RC25-010: Use clap conflicts_with for output args

**Epic**: REVIEW-CYCLE-2026-03-25
**Priority**: medium
**Depends on**: none
**Status**: done

## Goal

`--output` and `--output-dir` mutual exclusion is checked manually. Use clap's
`conflicts_with` attribute so `--help` shows the constraint and the framework
handles validation.

## Success Criteria

- [x] `conflicts_with` annotation is present on the relevant clap args.
- [x] Manual conflict check is removed.
- [x] `--help` shows the mutual exclusion.
- [x] All tests pass.
- [x] `/update` has been run after changes.
