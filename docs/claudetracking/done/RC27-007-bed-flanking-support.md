# RC27-007: Add BED flanking support

**Epic**: REVIEW-CYCLE-2026-03-27
**Priority**: medium
**Depends on**: none
**Status**: done

## Goal

The `--flank` flag currently applies to regions parsed from `--region` and
`--table` but not to regions parsed from `--bed`. Pass the flank value
through to `parse_regions_bed` so that BED regions are extended by the
requested amount on both sides, clamped to contig bounds.

## Success Criteria

- [x] `--bed` regions are extended by `--flank` bases on each side.
- [x] Flanked regions are clamped to `[0, contig_length]`.
- [x] A test verifies BED flanking produces the expected coordinates.
- [x] Existing non-BED flanking behaviour is unchanged.
- [x] All tests pass.
- [x] `/update` has been run after changes.
