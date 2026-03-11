# multiseqex

**MULTI SEQuence EXtractor** — a fast, parallel CLI tool for extracting multiple
sequences from FASTA files using `.fai` indexing.

Similar to `samtools faidx` but optimised for bulk extraction by leveraging
multiple CPU cores and flexible batch input formats (CSV/TSV tables with named
columns, list files, and inline regions).

## Usage

```text
Usage: multiseqex [OPTIONS] <FASTA>

Arguments:
  <FASTA>  Reference FASTA file (bgzipped ok if a matching .fai exists)

Options:
      --regions <REGIONS>        Comma-separated regions: chr:start-end, ...
      --list <LIST>              File with one region per line (chr:start-end)
      --table <TABLE>            CSV/TSV table with named columns (see below)
      --sv-table <SV_TABLE>      CSV/TSV SV table with named columns (see below)
      --flank <FLANK>            Flank size for position-mode tables
  -o, --output <OUTPUT>          Output FASTA file (default: stdout)
      --output-dir <OUTPUT_DIR>  Output directory (one file per region/SV pair)
      --threads <THREADS>        Number of worker threads (default: all CPUs)
      --no-build-fai             Error if .fai is missing instead of building it
  -h, --help                     Print help
  -V, --version                  Print version
```

### Examples

```bash
# Single region to stdout
multiseqex ref.fa --regions chr1:1000-2000

# Multiple regions to a file
multiseqex ref.fa --regions chr1:1000-2000,chr2:3000-4000 -o out.fa

# From a CSV table (range mode)
multiseqex ref.fa --table regions.csv -o out.fa

# From a CSV table (position mode with flanking)
multiseqex ref.fa --table positions.csv --flank 500 -o out.fa

# SV breakpoints to per-pair files
multiseqex ref.fa --sv-table variants.tsv --output-dir sv_seqs/

# One file per region
multiseqex ref.fa --table regions.csv --output-dir per_region/
```

## Table formats

Tables must have a **header row** with named columns. Column names are
**case-insensitive** and can appear in **any order**. Extra columns (e.g. `GENE`,
`STRAND`) are silently ignored.

### `--table`

| Column  | Required?                          | Description                      |
|---------|------------------------------------|----------------------------------|
| `CHROM` | Yes                                | Chromosome / contig name         |
| `START` | Yes (range mode)                   | 1-based inclusive start position |
| `END`   | Yes (range mode)                   | 1-based inclusive end position   |
| `POS`   | Yes (position mode, needs --flank) | Single coordinate position       |
| `NAME`  | No                                 | Region label for output naming   |

- **Range mode**: provide `CHROM`, `START`, `END`.
- **Position mode**: provide `CHROM`, `POS` and pass `--flank`.

### `--sv-table`

Each row produces **two regions** (left and right breakpoints).

| Column        | Required?                          | Description                |
|---------------|------------------------------------|----------------------------|
| `CHROM_LEFT`  | Yes                                | Left breakpoint chromosome |
| `START_LEFT`  | Yes (range mode)                   | Left breakpoint start      |
| `END_LEFT`    | Yes (range mode)                   | Left breakpoint end        |
| `POS_LEFT`    | Yes (position mode, needs --flank) | Left breakpoint position   |
| `CHROM_RIGHT` | Yes                                | Right breakpoint chromosome|
| `START_RIGHT` | Yes (range mode)                   | Right breakpoint start     |
| `END_RIGHT`   | Yes (range mode)                   | Right breakpoint end       |
| `POS_RIGHT`   | Yes (position mode, needs --flank) | Right breakpoint position  |
| `NAME`        | No                                 | SV identifier for naming   |

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

- Rust 1.85+ and Cargo (edition 2024)
- [samtools](http://www.htslib.org/doc/samtools.html) (optional — for
  pre-building `.fai` indexes)

If the FASTA file lacks a `.fai` index, `multiseqex` builds one automatically
(unless `--no-build-fai` is set).

## Documentation

See the [docs/](docs/) folder for detailed guides:

- [Usage guide](docs/usage.md) — full walkthrough of all input and output modes
- [Testing and benchmarking](docs/testing.md) — how to run tests and measure
  performance

## License

MIT
