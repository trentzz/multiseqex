# RC25-006: Add unit tests for core functions

**Epic**: REVIEW-CYCLE-2026-03-25
**Priority**: medium
**Depends on**: none
**Status**: done

## Goal

Zero unit tests exist. Add `#[cfg(test)]` module with unit tests for:
`parse_region_str`, `wrap_fasta`, `count_bases`, `build_fai`, `parse_fasta_header`,
and the table/region parsers. Cover expected paths, edge cases (empty input,
zero-length sequences, Windows line endings), and failure modes.

## Success Criteria

- [x] At least one unit test per function listed above.
- [x] Edge cases covered: empty input, zero-length sequence, malformed input.
- [x] All tests pass.
- [x] `/update` has been run after changes.
