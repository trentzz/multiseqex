# FEAT-002: Advanced Extraction Features

## Goal

Add features that make multiseqex substantially more useful for common
bioinformatics workflows.

## Motivation

Reverse complement extraction is essential for strand-aware analysis. Region
deduplication avoids wasted computation and confusing output. Sorted output
simplifies downstream processing. These are frequently requested capabilities
in FASTA extraction tools.

## Scope

- Reverse complement support (via a `--strand` or `STRAND` column).
- Region deduplication (optional, warn on duplicates).
- Sorted output by genomic coordinate (optional flag).

## Tasks

- F002-001: Add reverse complement support
- F002-002: Deduplicate identical regions
- F002-003: Add sorted output mode
