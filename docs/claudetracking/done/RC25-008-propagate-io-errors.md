# RC25-008: Propagate I/O errors in parse_regions_list

**Epic**: REVIEW-CYCLE-2026-03-25
**Priority**: medium
**Depends on**: none
**Status**: done

## Goal

`parse_regions_list` uses `.ok()?` inside `filter_map`, silently swallowing I/O
errors (e.g. invalid UTF-8). Propagate errors so users know when a line could
not be read.

## Success Criteria

- [x] I/O errors during region file parsing are reported to the user.
- [x] Silent line skipping on error no longer occurs.
- [x] All tests pass.
- [x] `/update` has been run after changes.
