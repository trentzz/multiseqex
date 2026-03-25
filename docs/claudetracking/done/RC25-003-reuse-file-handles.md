# RC25-003: Reuse file handles in parallel extraction

**Epic**: REVIEW-CYCLE-2026-03-25
**Priority**: high
**Depends on**: none
**Status**: done

## Goal

`extract_region` opens a new `File` for every region. With thousands of regions
and many threads, this exhausts file descriptors. Use a thread-local file handle
or open one handle per thread instead of one per region.

## Success Criteria

- [x] Each thread opens at most one file handle for the FASTA file.
- [x] Parallel extraction still works correctly.
- [x] All tests pass.
- [x] `/update` has been run after changes.
