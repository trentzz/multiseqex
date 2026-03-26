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

## Coordinate system

All coordinates are **1-based and inclusive on both ends**. A region
`chr1:100-200` extracts bases 100 through 200 inclusive (101 bases total).

When using position + flank syntax, `chr1:1000+500` means position 1000 with a
flank of 500 on each side. The resulting region is `chr1:500-1500` (1-based
inclusive). If the lower bound would fall below 1, it is clamped to 1. If the
upper bound exceeds the contig length, it is clamped to the contig length.

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

This extracts bases 500 to 1500 inclusive (position 1000 with 500 flanking bases each side).

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

An optional `NAME` column labels the region in FASTA headers (e.g. `>myregion
chr1:1000-2000`) and is used for filenames with `--output-dir` (e.g.
`myregion_1000_2000.fa`). Any additional columns (e.g. `GENE`, `STRAND`) are
silently ignored.
Column names are **case-insensitive** and can appear in any order.

### BED file (`--bed`)

A standard BED file (tab-separated). BED uses 0-based half-open coordinates.
They are converted to 1-based inclusive internally (start+1, end unchanged).
An optional fourth column provides a region name.

```bed
chr1	1000	2000
chr2	3000	4000	myregion
```

```bash
multiseqex ref.fa --bed regions.bed -o out.fa
```

Comments (lines starting with `#`) and blank lines are skipped.

The `--flank` flag extends each BED region symmetrically. Flanking is applied
after the coordinate conversion:

```bash
# Extend each region by 500bp on each side
multiseqex ref.fa --bed regions.bed --flank 500 -o out.fa
```

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

## Combining input sources

`--regions`, `--list`, and `--table` can be freely combined. All regions from
every source are merged into a single extraction. This is useful when you have
a base set of regions in a file but need to add a few extra on the command line.

```bash
# Combine a table with extra inline regions
multiseqex ref.fa --table regions.csv --regions chr5:100-200 -o out.fa

# Combine a list file with a table
multiseqex ref.fa --list regions.txt --table extra.csv -o out.fa

# All three at once
multiseqex ref.fa --regions chr1:1-500 --list regions.txt --table extra.csv -o out.fa
```

`--sv-table` can also be combined with the other sources. Regions from
`--sv-table` (left and right breakpoints) are merged alongside any regions from
`--regions`, `--list`, or `--table`.

## Output options

| Flag | Behaviour |
|------|-----------|
| *(none)* | Print FASTA to stdout |
| `-o out.fa` | Write all sequences to a single file |
| `--output-dir seqs/` | Write one `.fa` file per region (or per SV pair) |

`-o` and `--output-dir` cannot be used together.

When using `--output-dir` with `--sv-table`, each SV pair is written to one
file containing both breakpoint sequences.

## Reverse complement (`--rc`)

The `--rc` flag reverse-complements every extracted sequence before output.
All IUPAC ambiguity codes are supported.

```bash
multiseqex ref.fa --regions chr1:1000-2000 --rc -o out.fa
```

## Deduplication (`--dedup`)

The `--dedup` flag removes duplicate regions (same chromosome, start, end)
before extraction. The first occurrence of each region is kept. A message is
printed to stderr reporting how many duplicates were removed (unless `--quiet`
is set).

```bash
multiseqex ref.fa --table regions.csv --dedup -o out.fa
```

## Sorting (`--sort`)

The `--sort` flag sorts regions by chromosome (natural order: chr1, chr2, ...,
chr10) then by start position. This applies before extraction, so the output
follows sorted order.

```bash
multiseqex ref.fa --table regions.csv --sort -o out.fa
```

`--dedup` and `--sort` can be combined. Deduplication runs first.

## Quiet mode (`--quiet`)

The `-q` / `--quiet` flag suppresses all progress messages and warnings on
stderr. Error messages still appear.

```bash
multiseqex ref.fa --table big.csv -o out.fa --quiet
```

## Delimiter override (`--delimiter`)

By default, `--table` and `--sv-table` auto-detect the delimiter. Files with a
`.tsv` extension use tab. Other files are sniffed for tabs in the header line,
falling back to comma.

Use `--delimiter` to override:

```bash
multiseqex ref.fa --table data.txt --delimiter tab
multiseqex ref.fa --table data.txt --delimiter comma
multiseqex ref.fa --table data.txt --delimiter ";"
```

`--delimiter` requires `--table` or `--sv-table`.

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

## Performance notes

When writing to stdout or a single output file (`-o`), all extracted sequences
are buffered in memory to preserve input order. Memory usage therefore scales
with the total output size (number of regions multiplied by average region
length). For very large extraction jobs where memory is a concern, use
`--output-dir` instead. Each region is written independently in that mode,
keeping memory usage proportional to a single region at a time.
