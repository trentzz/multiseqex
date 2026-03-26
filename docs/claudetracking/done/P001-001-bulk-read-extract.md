# P001-001: Bulk-read optimisation for extract_region

**Epic**: PERF-001
**Priority**: medium
**Depends on**: none
**Status**: done

## Goal

Replace the line-by-line seek/read loop in `extract_region` with a single bulk
read. For a region spanning N FASTA lines, the current code does N seeks and N
reads. A single seek to the start byte, one read of the full byte range, and
in-memory newline stripping would reduce syscalls dramatically.

## Success Criteria

- [x] `extract_region` uses at most one seek and one read for any region.
- [x] Newline characters within the read buffer are stripped in-memory.
- [x] All existing tests pass (correctness unchanged).
- [x] Bulk-read groups coalesce nearby regions on the same contig (8 KB threshold).
- [x] /update has been run after changes.
