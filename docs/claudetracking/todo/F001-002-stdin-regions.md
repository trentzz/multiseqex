# F001-002: Accept region input from stdin

**Epic**: FEAT-001
**Priority**: low
**Depends on**: none
**Status**: todo

## Goal

Allow piping region input from other tools via stdin. When `--list -` is
passed, read regions from stdin instead of a file. This enables workflows
like `grep ... regions.bed | multiseqex ref.fa --list -`.

## Success Criteria

- [ ] `--list -` reads from stdin.
- [ ] Integration test pipes regions through stdin.
- [ ] Documentation updated.
- [ ] All existing tests pass.
- [ ] /update has been run after changes.
