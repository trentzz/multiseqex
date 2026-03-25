# RC25-001: Remove bogus bgzip claim

**Epic**: REVIEW-CYCLE-2026-03-25
**Priority**: high
**Depends on**: none
**Status**: done

## Goal

Remove the false claim that bgzipped FASTA files are supported. The code does
raw `File::open` + `Seek` + `read_exact` with no decompression. A bgzipped file
produces garbage silently. Remove the claim from CLI help (line 44) and README
(line 16). Optionally, detect bgzip magic bytes and exit with a clear error.

## Success Criteria

- [x] CLI help text no longer mentions bgzip.
- [x] README no longer mentions bgzip.
- [x] If a bgzipped file is passed, the tool either errors clearly or the docs
      state it is unsupported.
- [x] All tests pass.
- [x] `/update` has been run after changes.
