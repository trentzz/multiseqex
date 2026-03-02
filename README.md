# multiseqex
MULTI SEQuence EXtractor

`multiseqex` is a CLI tool that allows you to extract a sequence from a reference file similar to what `samtools faidx` does. The real advantage of this tool is being able to extract multiple sequences in a more efficient manner by leveraging multiple cores and batch processing.

## Usage

```text
Multi-sequence extractor for FASTA using FAI

Usage: multiseqex [OPTIONS] <FASTA>

Arguments:
  <FASTA>  Reference FASTA file (bgzipped FASTA is fine as long as a matching .fai exists)

Options:
      --regions <REGIONS>        Comma-separated regions: chr:start-end,chr2:start-end,...
      --list <LIST>              File with one region per line in the form chr:start-end
      --table <TABLE>            CSV/TSV with named columns (see Table Formats below)
      --sv-table <SV_TABLE>      CSV/TSV SV table with named columns (see Table Formats below)
      --flank <FLANK>            Flank size for position-mode tables (required with POS columns)
  -o, --output <OUTPUT>          Output FASTA file (single combined; default: stdout)
      --output-dir <OUTPUT_DIR>  Output directory for per-sequence FASTA files
      --threads <THREADS>        Number of worker threads (default: Rayon default, usually #cpus)
      --no-build-fai             Do not build a .fai if missing (error instead)
  -h, --help                     Print help
  -V, --version                  Print version
```

## Table Formats

Tables must have a **header row** with named columns. Column order does not matter, and extra columns (e.g. GENE, STRAND) are silently ignored. Column names are **case-insensitive**.

### `--table` (regular regions)

| Column  | Required?                         | Description                      |
|---------|-----------------------------------|----------------------------------|
| `CHROM` | Yes                               | Chromosome / contig name         |
| `START` | Yes (range mode)                  | 1-based inclusive start position |
| `END`   | Yes (range mode)                  | 1-based inclusive end position   |
| `POS`   | Yes (position mode, needs --flank)| Single coordinate position       |
| `NAME`  | No                                | Region label for output naming   |

**Mode detection:**
- **Range mode** (CHROM + START + END): extracts the given range directly.
- **Position mode** (CHROM + POS): requires `--flank`; extracts `POS - flank` to `POS + flank`.

Example TSV (range mode):
```
CHROM	START	END	GENE
chr1	1000	2000	TP53
chr2	5000	6000	BRCA1
```

Example CSV (position mode with `--flank 500`):
```
CHROM,POS,NAME
chr1,1500,region_A
chr2,5500,region_B
```

### `--sv-table` (structural variant regions)

| Column        | Required?                         | Description                |
|---------------|-----------------------------------|----------------------------|
| `CHROM_LEFT`  | Yes                               | Left breakpoint chromosome |
| `START_LEFT`  | Yes (range mode)                  | Left breakpoint start      |
| `END_LEFT`    | Yes (range mode)                  | Left breakpoint end        |
| `POS_LEFT`    | Yes (position mode, needs --flank)| Left breakpoint position   |
| `CHROM_RIGHT` | Yes                               | Right breakpoint chromosome|
| `START_RIGHT` | Yes (range mode)                  | Right breakpoint start     |
| `END_RIGHT`   | Yes (range mode)                  | Right breakpoint end       |
| `POS_RIGHT`   | Yes (position mode, needs --flank)| Right breakpoint position  |
| `NAME`        | No                                | SV identifier for naming   |

**Mode detection:**
- **Range mode** (all START/END columns present): extracts explicit ranges for both breakpoints.
- **Position mode** (POS_LEFT + POS_RIGHT): requires `--flank`; expands each position by the flank value.

Example TSV (range mode):
```
NAME	CHROM_LEFT	START_LEFT	END_LEFT	CHROM_RIGHT	START_RIGHT	END_RIGHT
SV001	chr1	1000	2000	chr3	5000	6000
SV002	chr2	3000	4000	chr5	7000	8000
```

Example CSV (position mode with `--flank 1000`):
```
CHROM_LEFT,POS_LEFT,CHROM_RIGHT,POS_RIGHT
chr1,1500,chr3,5500
```

## Installation

### From this github page

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

- Rust (version 1.60 or higher recommended)
- Cargo (Rust package manager, usually included with Rust)
- [samtools](http://www.htslib.org/doc/samtools.html) (optional, for building `.fai` index files)

If your FASTA file does not have a corresponding `.fai` index, `multiseqex` will attempt to build one automatically unless you use the `--no-build-fai` flag. Alternatively, you can manually generate the index using:

```bash
samtools faidx <your_fasta_file>
```
