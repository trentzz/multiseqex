# PERF-001: Extraction Performance

## Goal

Improve I/O performance for large-batch and large-region extractions.

## Motivation

The current `extract_region` function reads one line-width chunk at a time,
issuing a seek + read per FASTA line crossed. For a 1MB region with 60bp lines,
that is ~17,000 seek/read pairs. A bulk read with newline stripping would be
significantly faster.

Additionally, all extracted sequences are collected into memory before any
output is written. For large batches this wastes memory.

## Scope

- Bulk-read optimisation for `extract_region`.
- Streaming output to reduce peak memory.
- Benchmark before and after.

## Tasks

- P001-001: Bulk-read optimisation for extract_region
- P001-002: Streaming output to reduce memory usage
