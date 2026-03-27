# F001-002: Accept region input from stdin

**Epic**: FEAT-001
**Priority**: low
**Depends on**: none
**Status**: done

## Goal

Allow piping region input from other tools via stdin. When `--list -` is
passed, read regions from stdin instead of a file. This enables workflows
like `grep ... regions.bed | multiseqex ref.fa --list -`.

## Success Criteria

- [x] `--list -` reads from stdin.
- [x] Integration test pipes regions through stdin.
- [x] Documentation updated.
- [x] All existing tests pass.
- [x] /update has been run after changes.
