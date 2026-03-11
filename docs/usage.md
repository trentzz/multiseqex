# Usage Guide

`multiseqex` extracts one or more sequences from a FASTA file using its `.fai`
index. It supports four ways to specify regions, and can write output to stdout,
a single file, or a directory of per-region files.

## Quick start

```bash
# Extract a single region to stdout
multiseqex ref.fa --regions chr1:1000-2000

# Extract multiple regions to a file
multiseqex ref.fa --regions chr1:1000-2000,chr2:3000-4000 -o out.fa

# Extract from a CSV table
multiseqex ref.fa --table regions.csv -o out.fa
```

## Specifying regions

### Inline (`--regions`)

Comma-separated `chr:start-end` strings:

```bash
multiseqex ref.fa --regions chr1:100-200,chr2:300-400
```

Position + flank syntax is also supported:

```bash
multiseqex ref.fa --regions chr1:1000+500
```

This extracts bases 500–1500 (position 1000 ± 500).

### List file (`--list`)

A plain-text file with one region per line. Blank lines and lines starting with
`#` are ignored.

```
# my regions
chr1:100-200
chr2:300-400
chr3:500-600
```

```bash
multiseqex ref.fa --list regions.txt
```

### CSV/TSV table (`--table`)

A delimited file with **named column headers**. The delimiter is auto-detected
from the file extension (`.tsv` → tab, anything else → comma).

**Range mode** — requires `CHROM`, `START`, `END`:

```csv
CHROM,START,END
chr1,1000,2000
chr2,3000,4000
```

**Position mode** — requires `CHROM`, `POS`, and the `--flank` flag:

```csv
CHROM,POS,NAME
chr1,1500,regionA
chr2,3500,regionB
```

```bash
multiseqex ref.fa --table positions.csv --flank 500
```

An optional `NAME` column labels the region for output file naming.
Any additional columns (e.g. `GENE`, `STRAND`) are silently ignored.
Column names are **case-insensitive** and can appear in any order.

### SV table (`--sv-table`)

For structural variants, each row produces **two regions** (left and right
breakpoints).

**Range mode** — requires `CHROM_LEFT`, `START_LEFT`, `END_LEFT`, `CHROM_RIGHT`,
`START_RIGHT`, `END_RIGHT`:

```tsv
NAME	CHROM_LEFT	START_LEFT	END_LEFT	CHROM_RIGHT	START_RIGHT	END_RIGHT
SV001	chr1	1000	2000	chr3	5000	6000
```

**Position mode** — requires `CHROM_LEFT`, `POS_LEFT`, `CHROM_RIGHT`,
`POS_RIGHT`, plus `--flank`:

```csv
CHROM_LEFT,POS_LEFT,CHROM_RIGHT,POS_RIGHT
chr1,1500,chr3,5500
```

```bash
multiseqex ref.fa --sv-table sv_positions.csv --flank 1000
```

## Output options

| Flag | Behaviour |
|------|-----------|
| *(none)* | Print FASTA to stdout |
| `-o out.fa` | Write all sequences to a single file |
| `--output-dir seqs/` | Write one `.fa` file per region (or per SV pair) |

`-o` and `--output-dir` cannot be used together.

When using `--output-dir` with `--sv-table`, each SV pair is written to one
file containing both breakpoint sequences.

## FAI index

`multiseqex` requires a `.fai` index alongside the FASTA file. If one is not
found, it is built automatically. To suppress auto-building (e.g. for read-only
filesystems), use `--no-build-fai`.

You can also pre-build the index with samtools:

```bash
samtools faidx ref.fa
```

## Threading

By default all available CPU cores are used. Override with `--threads`:

```bash
multiseqex ref.fa --table big_table.csv --threads 4
```
