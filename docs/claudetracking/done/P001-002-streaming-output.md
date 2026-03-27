# P001-002: Streaming output to reduce memory usage

**Epic**: PERF-001
**Priority**: low
**Depends on**: P001-001
**Status**: done

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

- [x] `--output-dir` mode writes files as they are extracted (no full collect).
- [x] Single-file/stdout mode uses ordered-slot streaming with bulk-read groups.
- [x] Output order matches input region order (deterministic).
- [x] All existing tests pass.
- [x] /update has been run after changes.
