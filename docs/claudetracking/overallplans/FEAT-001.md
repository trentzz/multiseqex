# FEAT-001: Input Format and Usability Features

## Goal

Expand input format support and improve usability for common bioinformatics
workflows.

## Motivation

BED is the most widely used interval format in genomics. Not supporting it is a
significant gap. Stdin support for region input enables piping from other tools.
Progress reporting helps users working with large region sets.

## Scope

- BED file support (0-based, half-open).
- Stdin support for region input.
- Progress indicator for large extractions.
- Validate 1-based coordinates (reject 0).

## Tasks

- F001-001: Add BED file input support
- F001-002: Accept region input from stdin
- F001-003: Add progress reporting for large extractions
- F001-004: Reject zero-start coordinates in 1-based mode
