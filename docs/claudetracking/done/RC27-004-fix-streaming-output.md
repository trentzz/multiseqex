# RC27-004: Fix streaming output to be truly streaming

**Epic**: REVIEW-CYCLE-2026-03-27
**Priority**: high
**Depends on**: none
**Status**: done

## Goal

The current output path buffers all extracted sequences in memory before
writing. For large extraction jobs this causes unnecessary memory pressure.
Refactor the output to write each sequence as it is extracted, streaming
results to stdout or the output file.

If truly streaming output is too complex to implement safely (e.g. due to
parallel extraction ordering requirements), at minimum document the
limitation clearly in the CLI help and usage docs so users are aware of
memory requirements for large jobs.

## Success Criteria

- [x] Output is written incrementally as sequences are extracted, or the
      limitation is documented in CLI help and usage docs.
- [x] Memory usage does not grow proportionally to the number of regions
      (verified by inspection or a simple test).
- [x] Existing output is byte-identical before and after the change.
- [x] All tests pass.
- [x] `/update` has been run after changes.
