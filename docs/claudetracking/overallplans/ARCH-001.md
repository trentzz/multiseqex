# ARCH-001: Module Decomposition

## Goal

Split `src/main.rs` (1124 lines) into logical modules. The single-file
architecture makes the codebase harder to navigate and test in isolation.

## Motivation

All functionality (CLI parsing, FAI building/reading, region parsing, sequence
extraction, output writing) lives in one file. Separating concerns into modules
improves readability, enables targeted unit testing, and prepares the codebase
for library extraction.

## Scope

- Extract FAI logic into `src/fai.rs`.
- Extract region parsing into `src/regions.rs`.
- Extract sequence extraction into `src/extract.rs`.
- Extract output writing into `src/output.rs`.
- Keep CLI definition and `main()` in `src/main.rs`.
- All existing tests must continue to pass.

## Tasks

- A001-001: Split main.rs into modules
