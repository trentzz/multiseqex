# F001-005: Delimiter detection by content sniffing

**Epic**: FEAT-001
**Priority**: medium
**Depends on**: none
**Status**: done

## Goal

Currently delimiter detection is extension-only (.tsv = tab, else comma). Sniff
the header line for tabs as a fallback when the extension is ambiguous. Add a
`--delimiter` CLI flag to override detection entirely.

## Success Criteria

- [x] If a file has no .tsv extension, the first line is checked for tabs.
- [x] `--delimiter` flag accepts "tab", "comma", or a literal character.
- [x] Existing table tests still pass.
- [x] New tests cover: .txt file with tabs detected correctly, --delimiter override.
- [x] `/update` has been run after changes.
