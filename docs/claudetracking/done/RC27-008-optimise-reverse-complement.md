# RC27-008: Optimise reverse_complement to byte-level iteration

**Epic**: REVIEW-CYCLE-2026-03-27
**Priority**: medium
**Depends on**: none
**Status**: done

## Goal

The current `reverse_complement` function works at the `char` level, which
involves unnecessary UTF-8 decoding for what is purely ASCII data. Refactor
to iterate over bytes directly, complement using a lookup table or match on
`u8`, and reverse the byte buffer. This avoids per-character decode overhead
on large sequences.

## Success Criteria

- [x] `reverse_complement` operates on bytes, not chars.
- [x] Output is identical to the previous implementation for all valid bases.
- [x] Invalid bytes are handled the same way as before (or better).
- [x] Unit tests confirm correctness for standard and ambiguous IUPAC bases.
- [x] All tests pass.
- [x] `/update` has been run after changes.
