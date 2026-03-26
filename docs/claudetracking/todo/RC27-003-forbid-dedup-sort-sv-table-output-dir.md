# RC27-003: Forbid --dedup/--sort with --sv-table --output-dir

**Epic**: REVIEW-CYCLE-2026-03-27
**Priority**: high
**Depends on**: none
**Status**: todo

## Goal

When `--sv-table` is used with `--output-dir`, each sample produces a
separate output file with paired regions. The `--dedup` and `--sort` flags
operate on a flat list of regions and break the SV pairing semantics (pairs
may be deduplicated away or reordered so forward/reverse no longer
correspond).

Either forbid the combination with a clear error, or implement
pairing-aware dedup and sort. Forbidding is the safer initial approach.

## Success Criteria

- [ ] `--dedup` or `--sort` with `--sv-table --output-dir` produces a clear
      error message.
- [ ] The combination still works if `--output-dir` is not used (flat output
      mode), or is also forbidden with rationale documented.
- [ ] Tests cover the rejected combinations.
- [ ] All tests pass.
- [ ] `/update` has been run after changes.
