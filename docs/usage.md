# Usage Guide

`multiseqex` extracts one or more sequences from FASTA files using `.fai`
indexing. It supports many input formats, sequence transforms, masking, interval
operations, and multiple output modes.

## Quick start

```bash
# Extract a single region to stdout
multiseqex ref.fa --regions chr1:1000-2000

# Extract from a BED file
multiseqex ref.fa --bed regions.bed -o out.fa

# VCF variant context with flanking
multiseqex ref.fa --vcf variants.vcf --flank 100 -o out.fa

# GFF gene extraction
multiseqex ref.fa --gff annotations.gff3 -o genes.fa

# Per-region statistics
multiseqex ref.fa --bed regions.bed --stats
```

## Coordinate systems

All internal coordinates are **1-based inclusive**. A region `chr1:100-200`
extracts bases 100 through 200 inclusive (101 bases).

**BED files** use 0-based half-open coordinates. `multiseqex` converts them
automatically: BED `chr1 999 2000` becomes the internal region `chr1:1000-2000`.

**Position + flank** syntax (`chr1:1000+500`) expands to `chr1:500-1500`. If the
lower bound falls below 1, it is clamped to 1. If the upper bound exceeds the
contig length, it is clamped to the contig length.

## Input formats

### Inline regions (`--regions`)

Comma-separated `chr:start-end` or `chr:pos+flank` strings:

```bash
multiseqex ref.fa --regions chr1:100-200,chr2:300-400
multiseqex ref.fa --regions chr1:1000+500
```

### List file (`--list`)

A plain-text file with one region per line. Blank lines and lines starting with
`#` are ignored. Use `-` to read from stdin.

```
# my regions
chr1:100-200
chr2:300-400
```

```bash
multiseqex ref.fa --list regions.txt
cat regions.txt | multiseqex ref.fa --list -
```

### BED file (`--bed`)

Standard BED format (tab-separated). Coordinates are 0-based half-open. An
optional fourth column provides a region name. Column 6 (strand) is recognised
when present.

```bed
chr1	999	2000	regionA
chr2	2999	4000	regionB	0	+
```

```bash
multiseqex ref.fa --bed regions.bed -o out.fa
```

Comments (`#`) and blank lines are skipped. Empty intervals (start == end) and
malformed coordinates (start > end) are rejected.

### CSV/TSV table (`--table`)

A delimited file with a header row. Column names are case-insensitive and can
appear in any order. Extra columns are ignored.

The delimiter is auto-detected: `.tsv` files use tab; other files are sniffed
for tabs, falling back to comma. Use `--delimiter` to override.

**Range mode** requires `CHROM`, `START`, `END`:

```csv
CHROM,START,END,NAME
chr1,1000,2000,regionA
chr2,3000,4000,regionB
```

**Position mode** requires `CHROM`, `POS`, and the `--flank` flag:

```csv
CHROM,POS,NAME
chr1,1500,regionA
chr2,3500,regionB
```

```bash
multiseqex ref.fa --table positions.csv --flank 500
```

Optional columns: `NAME` (region label) and `STRAND` (`+`, `-`, or `.`).

### SV table (`--sv-table`)

Each row produces two regions (left and right breakpoints).

**Range mode** requires `CHROM_LEFT`, `START_LEFT`, `END_LEFT`, `CHROM_RIGHT`,
`START_RIGHT`, `END_RIGHT`:

```tsv
NAME	CHROM_LEFT	START_LEFT	END_LEFT	CHROM_RIGHT	START_RIGHT	END_RIGHT
SV001	chr1	1000	2000	chr3	5000	6000
```

**Position mode** requires `CHROM_LEFT`, `POS_LEFT`, `CHROM_RIGHT`, `POS_RIGHT`,
plus `--flank`:

```bash
multiseqex ref.fa --sv-table sv_positions.csv --flank 1000
```

With `--output-dir`, each SV pair is written to a single file containing both
breakpoint sequences.

### VCF (`--vcf`)

Extracts a region spanning `POS` to `POS + len(REF) - 1` for each VCF record.
The `ID` field is used as the region name (unless `.`). `REF` and `ALT` are
included in the FASTA header description.

```bash
multiseqex ref.fa --vcf variants.vcf -o out.fa
multiseqex ref.fa --vcf variants.vcf --flank 100 -o context.fa
```

Flanking extends each variant region symmetrically (or asymmetrically with
`--flank-left` / `--flank-right`).

### GFF3/GTF (`--gff`)

Extracts regions for features matching `--gff-feature` (default: `gene`).

```bash
multiseqex ref.fa --gff annotations.gff3 -o genes.fa
multiseqex ref.fa --gff annotations.gtf --gff-feature exon -o exons.fa
```

Flanking can be applied to GFF regions using `--flank`, `--flank-left`, or
`--flank-right`.

### Whole contigs (`--contigs`, `--contig-list`)

Extract entire contigs by name.

```bash
# Inline, comma-separated
multiseqex ref.fa --contigs chr1,chr2,chrX

# From a file (one name per line; # comments and blank lines skipped)
multiseqex ref.fa --contig-list contigs.txt
```

### Multiple FASTA files

Pass multiple FASTA files as positional arguments. Contigs are looked up across
all files. Each contig must appear in exactly one file.

```bash
multiseqex genome_part1.fa genome_part2.fa --regions chr1:1000-2000,chr5:500-600
```

## Region manipulation

### Flanking (`--flank`, `--flank-left`, `--flank-right`)

`--flank N` extends regions by N bases on each side. For asymmetric flanking,
use `--flank-left` and `--flank-right` together.

```bash
# Symmetric: 500bp each side
multiseqex ref.fa --bed regions.bed --flank 500

# Asymmetric: 200bp left, 800bp right
multiseqex ref.fa --bed regions.bed --flank-left 200 --flank-right 800
```

Flanking applies to BED, table (position mode), VCF, and GFF inputs. Bounds are
clamped to contig boundaries.

### Deduplication (`--dedup`)

Removes duplicate regions (same chromosome, start, end) before extraction. The
first occurrence is kept. Reports the count of removed duplicates on stderr.

```bash
multiseqex ref.fa --table regions.csv --dedup -o out.fa
```

### Sorting (`--sort`)

Sorts regions by natural chromosome order (chr1, chr2, ..., chr10, ...) then by
start position.

```bash
multiseqex ref.fa --table regions.csv --sort -o out.fa
```

`--dedup` and `--sort` can be combined. Deduplication runs first.

### Merging (`--merge`, `--merge-distance`)

Merges overlapping or book-ended regions on the same chromosome. Implies
`--sort`. Use `--merge-distance` to merge regions within a given gap (default:
0).

```bash
multiseqex ref.fa --bed regions.bed --merge -o out.fa
multiseqex ref.fa --bed regions.bed --merge --merge-distance 100 -o out.fa
```

### Interval operations (`--subtract`, `--intersect`)

`--intersect` keeps only the portions of input regions that overlap with a BED
file. `--subtract` removes the overlapping portions, potentially splitting
regions.

Intersection runs before subtraction when both are given.

```bash
# Keep only regions overlapping with targets.bed
multiseqex ref.fa --bed regions.bed --intersect targets.bed -o out.fa

# Remove repeats from regions
multiseqex ref.fa --bed regions.bed --subtract repeats.bed -o out.fa
```

### K-mer tiling (`--tile`, `--step`)

Tiles each region into fixed-width windows. `--step` sets the stride (default:
same as `--tile`, producing non-overlapping tiles). The final tile in each
region may be shorter than `--tile`.

```bash
# Non-overlapping 100bp tiles
multiseqex ref.fa --bed regions.bed --tile 100 -o tiles.fa

# Overlapping 100bp tiles with 50bp step
multiseqex ref.fa --bed regions.bed --tile 100 --step 50 -o tiles.fa
```

## Output formats

### FASTA (default)

Standard FASTA output. Line width defaults to 60 characters.

```bash
# Custom line width
multiseqex ref.fa --regions chr1:1-1000 --line-width 80

# No wrapping (single line per sequence)
multiseqex ref.fa --regions chr1:1-1000 --no-wrap
```

### FASTQ (`--fastq`, `--qual`)

Produces FASTQ output with a constant quality character (default: `I`, phred
40).

```bash
multiseqex ref.fa --bed regions.bed --fastq -o out.fq
multiseqex ref.fa --bed regions.bed --fastq --qual "F" -o out.fq
```

### TSV (`--tab-out`)

Tab-separated output with columns: chr, start, end, name, sequence.

```bash
multiseqex ref.fa --bed regions.bed --tab-out -o out.tsv
```

### Statistics (`--stats`)

Prints a TSV table with per-region statistics instead of extracting sequences.
Columns: chr, start, end, name, length, gc_percent, n_count, masked_count.

```bash
multiseqex ref.fa --bed regions.bed --stats
multiseqex ref.fa --bed regions.bed --stats > stats.tsv
```

### Output destinations

| Flag | Behaviour |
|------|-----------|
| *(none)* | Print to stdout |
| `-o out.fa` | Write all sequences to a single file |
| `--output-dir seqs/` | Write one file per region (or per SV pair) |

`-o` and `--output-dir` cannot be used together.

### Custom headers (`--name-template`)

Format FASTA/FASTQ headers using placeholders:

```bash
multiseqex ref.fa --bed regions.bed --name-template "{chr}_{start}_{end}" -o out.fa
```

Available placeholders: `{chr}`, `{start}`, `{end}`, `{name}`, `{length}`,
`{index}`, `{strand}`.

## Sequence transforms

### Reverse complement (`--rc`)

Reverse-complements every extracted sequence. All IUPAC ambiguity codes are
supported. Strand-aware: BED column 6 strand is applied automatically, and
`--rc` is XOR'd with the strand annotation.

```bash
multiseqex ref.fa --regions chr1:1000-2000 --rc -o out.fa
```

### DNA to RNA (`--to-rna`)

Converts T to U in output sequences.

```bash
multiseqex ref.fa --regions chr1:1000-2000 --to-rna
```

### Translation (`--translate`)

Translates to amino acids using the standard genetic code. Reading frame starts
at position 1. Stop codons are represented as `*`.

```bash
multiseqex ref.fa --regions chr1:1000-2000 --translate
```

### Case conversion (`--uppercase`, `--lowercase`)

Forces all output bases to uppercase or lowercase.

```bash
multiseqex ref.fa --bed regions.bed --uppercase -o out.fa
```

## Masking

Mask bases within extracted sequences using a BED file of mask regions.

```bash
# Hard mask (replace with N, the default)
multiseqex ref.fa --bed regions.bed --mask-bed repeats.bed -o out.fa

# Soft mask (lowercase)
multiseqex ref.fa --bed regions.bed --mask-bed repeats.bed --soft-mask -o out.fa
```

`--hard-mask` is the default when `--mask-bed` is given. Use `--soft-mask` for
lowercase masking instead.

## Alternate allele sequences (`--alt-seq`, `--alt-seq-both`)

These flags generate sequences where the REF allele is replaced with the ALT
allele, producing variant-modified output.

### With `--vcf`

Each VCF record's reference context is extracted (optionally with `--flank`),
then the REF bases are replaced by the ALT allele:

```bash
# SNP: REF=C ALT=G at chr1:5 with 3bp flanking
# Reference context: AACCCGG -> Alternate: AACGCGG
multiseqex ref.fa --vcf variants.vcf --flank 3 --alt-seq -o alt.fa
```

### With `--table`

The table must contain `REF` and `ALT` columns alongside `CHROM` and `POS`.
Position mode with `--flank` is required:

```bash
multiseqex ref.fa --table variants.csv --flank 3 --alt-seq -o alt.fa
```

### Both reference and alternate (`--alt-seq-both`)

`--alt-seq-both` implies `--alt-seq`. For each variant, two entries are emitted:
the unmodified reference sequence (tagged `ref_seq` in the header) followed by
the alternate sequence (tagged `alt_seq`):

```bash
multiseqex ref.fa --vcf variants.vcf --flank 100 --alt-seq-both -o both.fa
```

### Multi-allelic sites

When a VCF record or table row has a comma-separated ALT field (e.g. `G,C`),
each alternative allele produces a separate output entry:

```bash
# REF=T ALT=G,C -> two output sequences, one for each ALT
multiseqex ref.fa --vcf multi.vcf --flank 3 --alt-seq -o multi_alt.fa
```

### Variant types

- **SNP**: single-base replacement. Output length equals the reference context
  length.
- **Insertion**: ALT is longer than REF. Output is longer than the reference
  context.
- **Deletion**: ALT is shorter than REF. Output is shorter than the reference
  context.

### Constraints

- `--alt-seq` requires `--vcf` or `--table` (with REF/ALT columns). Using it
  with `--regions`, `--bed`, or other input sources alone produces an error.
- `--alt-seq` and `--alt-seq-both` conflict with `--sv-table`.

## Combining input sources

Most input flags can be combined freely. All regions from every source are
collected and extracted together.

```bash
# Table + inline regions
multiseqex ref.fa --table regions.csv --regions chr5:100-200 -o out.fa

# List + BED + VCF
multiseqex ref.fa --list regions.txt --bed extra.bed --vcf variants.vcf -o out.fa
```

`--sv-table` conflicts with `--regions`, `--table`, `--list`, `--contigs`,
`--contig-list`, `--vcf`, and `--gff`. It must be used alone or combined only
with `--bed`.

## Indexing and streaming

### FAI index

`multiseqex` requires a `.fai` index alongside each FASTA file. If one is not
found, it is built automatically. Use `--no-build-fai` to error instead (useful
for read-only filesystems).

Pre-build with samtools:

```bash
samtools faidx ref.fa
```

### No-index mode (`--no-index`)

Loads the entire FASTA into memory without an FAI index. Required when reading
from stdin. Supports only a single FASTA file.

```bash
cat ref.fa | multiseqex - --no-index --regions chr1:1000-2000
```

### Bgzip support

Bgzipped and gzipped FASTA files are decompressed transparently. The `.fai`
index is looked up next to the original compressed file first, then next to the
decompressed temporary file.

```bash
multiseqex ref.fa.gz --regions chr1:1000-2000
```

## Performance notes

- **Threading**: all available CPU cores are used by default. Override with
  `--threads N`.
- **Bulk-read optimisation**: nearby regions on the same contig are read in a
  single I/O operation, reducing seek overhead.
- **Streaming output**: stdout and single-file output (`-o`) buffer results in
  memory to preserve input order while extracting in parallel. Memory usage
  scales with total output size.
- **Per-region output**: `--output-dir` writes each region independently,
  keeping memory usage proportional to a single region at a time.
- **Progress bar**: shown on stderr when writing to a file (suppressed by
  `--quiet`).

## Limitations

- `--stats` is not supported with `--no-index`.
- `--sv-table` with `--output-dir` cannot be combined with `--merge`, `--dedup`,
  `--sort`, or `--tile` because these break the paired region invariant.
- `--no-index` supports only a single FASTA file (or stdin via `-`).
- `--to-rna` and `--translate` conflict with `--stats`.
- `--flank-left` and `--flank-right` must be specified together.

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

## Quiet mode (`--quiet`)

The `-q` / `--quiet` flag suppresses all progress messages, warnings, and the
progress bar on stderr. Error messages still appear.

```bash
multiseqex ref.fa --table big.csv -o out.fa --quiet
```
