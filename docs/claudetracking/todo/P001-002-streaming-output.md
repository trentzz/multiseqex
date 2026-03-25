# P001-002: Streaming output to reduce memory usage

**Epic**: PERF-001
**Priority**: low
**Depends on**: P001-001
**Status**: todo

## Goal

Currently all extracted sequences are collected into a `Vec<(Region, String)>`
before any output is written. For very large extractions (e.g. 100k regions of
1MB each), this holds all sequence data in memory simultaneously.

For single-file and stdout output, sequences must be written in input order
(parallel extraction can reorder). Investigate whether an ordered parallel
write is feasible, or whether a two-pass approach (parallel extract to temp
files, sequential merge) is better.

For `--output-dir` mode, each file is independent and could be written
immediately after extraction.

## Success Criteria

- [ ] `--output-dir` mode writes files as they are extracted (no full collect).
- [ ] Single-file mode either streams or documents the memory constraint.
- [ ] All existing tests pass.
- [ ] /update has been run after changes.
