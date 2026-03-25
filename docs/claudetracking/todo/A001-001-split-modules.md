# A001-001: Split main.rs into modules

**Epic**: ARCH-001
**Priority**: medium
**Depends on**: none
**Status**: todo

## Goal

Decompose the 1124-line `src/main.rs` into focused modules:

- `src/fai.rs`: `FaiRecord`, `fai_path_for`, `build_fai`, `read_fai`.
- `src/regions.rs`: `Region`, all `parse_regions_*` functions, `TableMode`, `SvMode`.
- `src/extract.rs`: `extract_region`.
- `src/output.rs`: `wrap_fasta`, `write_sequences`, `write_per_file`, `write_sv_per_file`.
- `src/main.rs`: `Cli`, `main()`, `detect_gzip_and_reject`, `validate_and_clamp_regions`.

Move unit tests alongside their modules. Integration tests remain in
`tests/integration.rs`.

## Success Criteria

- [ ] `src/main.rs` is under 150 lines.
- [ ] Each module has its own file with relevant unit tests.
- [ ] All 33 existing tests pass without modification.
- [ ] `cargo clippy -- -D warnings` is clean.
- [ ] /update has been run after changes.
