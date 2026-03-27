# multiseqex

**MULTI SEQuence EXtractor** — a fast, parallel CLI tool for extracting
sequences from FASTA files using `.fai` indexing.

Similar to `samtools faidx` but built for bulk extraction. Supports multiple
input formats (BED, VCF, GFF, CSV/TSV tables, inline regions), sequence
transforms (reverse complement, RNA conversion, translation), masking,
interval arithmetic, and parallel output across multiple CPU cores.

## Installation

### From crates.io

```bash
cargo install multiseqex
```

### From GitHub

```bash
cargo install --git https://github.com/trentzz/multiseqex
```

### Build from source

```bash
git clone https://github.com/trentzz/multiseqex.git
cd multiseqex
cargo build --release
cp target/release/multiseqex ~/.local/bin/
```

### Prerequisites

- Rust 1.87+ (edition 2024)
- [samtools](http://www.htslib.org/doc/samtools.html) (optional, for
  pre-building `.fai` indexes)

If the FASTA file lacks a `.fai` index, `multiseqex` builds one automatically
(unless `--no-build-fai` is set).

## Quick start

```bash
# Single region to stdout
multiseqex ref.fa --regions chr1:1000-2000

# Multiple regions from a BED file
multiseqex ref.fa --bed regions.bed -o out.fa

# VCF variant context extraction
multiseqex ref.fa --vcf variants.vcf --flank 100 -o out.fa

# GFF gene extraction
multiseqex ref.fa --gff annotations.gff3 --gff-feature gene -o genes.fa

# VCF alternate allele sequences with flanking context
multiseqex ref.fa --vcf variants.vcf --flank 100 --alt-seq -o alt.fa

# Both reference and alternate sequences
multiseqex ref.fa --vcf variants.vcf --flank 100 --alt-seq-both -o both.fa

# Region statistics (GC%, length, masking)
multiseqex ref.fa --bed regions.bed --stats

# Translate extracted sequences to amino acids
multiseqex ref.fa --regions chr1:1000-2000 --translate

# K-mer tiling with 100bp windows, 50bp step
multiseqex ref.fa --bed regions.bed --tile 100 --step 50 -o tiles.fa
```

## Options reference

### Input formats

| Flag | Description |
|------|-------------|
| `<FASTA>...` | One or more reference FASTA files (positional). Supports bgzip/gzip (transparent decompression). Use `-` for stdin with `--no-index`. |
| `--regions` | Comma-separated regions: `chr:start-end`, `chr:pos+flank` |
| `--list` | File with one region per line. Use `-` for stdin. |
| `--bed` | BED file (0-based half-open). Supports optional name (col 4) and strand (col 6). |
| `--table` | CSV/TSV with header. Columns: CHROM, START, END (range) or CHROM, POS (position + `--flank`). Optional: NAME, STRAND. |
| `--sv-table` | SV paired-region table. Columns: CHROM_LEFT/RIGHT, START/END_LEFT/RIGHT or POS_LEFT/RIGHT. |
| `--vcf` | VCF file. Extracts REF span per record. ID used as name; REF/ALT in header. |
| `--gff` | GFF3/GTF annotation file. Use `--gff-feature` to filter (default: `gene`). |
| `--alt-seq` | Generate alternate-allele sequences. Requires `--vcf` or `--table` with REF/ALT columns. |
| `--alt-seq-both` | Output both reference and alternate sequences. Implies `--alt-seq`. |
| `--contigs` | Comma-separated contig names to extract in full. |
| `--contig-list` | File with one contig name per line to extract in full. |

### Region manipulation

| Flag | Description |
|------|-------------|
| `--flank` | Symmetric flank size for position-mode regions. |
| `--flank-left` | Left-side flank (must pair with `--flank-right`). |
| `--flank-right` | Right-side flank (must pair with `--flank-left`). |
| `--dedup` | Remove duplicate regions (same chr, start, end). |
| `--sort` | Sort by natural chromosome order then start position. |
| `--merge` | Merge overlapping/book-ended regions. Implies `--sort`. |
| `--merge-distance` | Maximum gap for merging (default: 0). Requires `--merge`. |
| `--subtract` | BED file of intervals to subtract from input regions. |
| `--intersect` | BED file of intervals to intersect with input regions. |
| `--tile` | Tile each region into windows of this size (bases). |
| `--step` | Step size for tiling (default: same as `--tile`). |

### Output

| Flag | Description |
|------|-------------|
| `-o, --output` | Write all sequences to a single file (default: stdout). |
| `--output-dir` | Write one file per region (or per SV pair). |
| `--line-width` | FASTA line width (default: 60). Set to 0 to disable wrapping. |
| `--no-wrap` | Disable FASTA line wrapping (shorthand for `--line-width 0`). |
| `--tab-out` | TSV output: chr, start, end, name, sequence. |
| `--fastq` | FASTQ output with constant quality character. |
| `--qual` | Quality character for FASTQ (default: `I`, phred 40). |
| `--stats` | Print per-region statistics (TSV) instead of sequences. |
| `--name-template` | Custom header format. Placeholders: `{chr}`, `{start}`, `{end}`, `{name}`, `{length}`, `{index}`, `{strand}`. |
| `--rc` | Reverse complement all extracted sequences. |

### Transforms

| Flag | Description |
|------|-------------|
| `--to-rna` | Convert T to U (DNA to RNA). |
| `--translate` | Translate to amino acids (standard genetic code). Stop codons as `*`. |
| `--uppercase` | Force all output bases to uppercase. |
| `--lowercase` | Force all output bases to lowercase. |

### Masking

| Flag | Description |
|------|-------------|
| `--mask-bed` | BED file defining regions to mask within extracted sequences. |
| `--hard-mask` | Replace masked bases with N (default when `--mask-bed` is given). |
| `--soft-mask` | Lowercase masked bases instead of replacing with N. |

### Other

| Flag | Description |
|------|-------------|
| `--delimiter` | Override delimiter for `--table` / `--sv-table`. Accepts `tab`, `comma`, or a single character. |
| `--threads` | Number of worker threads (default: all available CPUs). |
| `--no-build-fai` | Error if `.fai` is missing instead of building one. |
| `--no-index` | Scan FASTA sequentially without an FAI index. Loads into memory. Required for stdin (`-`). |
| `-q, --quiet` | Suppress progress messages, warnings, and the progress bar. |

## Feature highlights

- **Multiple FASTA files**: pass several FASTA files as positional arguments.
  Contigs are looked up across all files (each contig must appear in exactly
  one file).
- **Bgzip support**: bgzipped and gzipped FASTA files are decompressed
  transparently.
- **Progress bar**: shown on stderr when writing to a file (unless `--quiet`).
- **Bulk-read optimisation**: nearby regions on the same contig are read in a
  single I/O operation, reducing seek overhead.
- **Streaming output**: stdout and single-file output buffer results in memory
  to preserve input order while extracting in parallel.
- **Interval arithmetic**: `--subtract` and `--intersect` apply set operations
  against a BED file before extraction.
- **K-mer tiling**: `--tile` and `--step` break regions into fixed-width
  windows for downstream analysis.
- **Coordinate systems**: inline regions and tables use 1-based inclusive
  coordinates. BED uses 0-based half-open (converted internally).

## Documentation

See the [docs/](docs/) folder for detailed guides:

- [Usage guide](docs/usage.md) — full walkthrough of all input and output modes
- [Testing and benchmarking](docs/testing.md) — how to run tests and measure
  performance

## Version

Current release: **0.2.1**
MSRV: **1.87** (Rust edition 2024)

## Licence

MIT
