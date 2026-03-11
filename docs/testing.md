# Testing and Benchmarking

## Running tests

The project uses Rust's built-in test framework with
[`assert_cmd`](https://crates.io/crates/assert_cmd) for integration testing.

```bash
# Run all tests
cargo test

# Run with output shown
cargo test -- --nocapture

# Run a specific test
cargo test table_csv_range
```

## Test structure

```
tests/
├── integration.rs          # Integration tests (CLI-level)
└── fixtures/
    ├── test.fa             # Small FASTA (3 contigs: chr1=204bp, chr2=140bp, chr3=72bp)
    ├── test.fa.fai         # Corresponding FAI index
    ├── regions.txt         # List-file input
    ├── table_range.csv     # CSV table (range mode)
    ├── table_range.tsv     # TSV table (range mode, with extra GENE column)
    ├── table_pos.csv       # CSV table (position mode)
    ├── sv_table_range.tsv  # SV table (range mode)
    └── sv_table_pos.csv    # SV table (position mode)
```

## What the tests cover

| Category | Tests |
|----------|-------|
| `--regions` | Single region, multiple regions, clamping to contig length |
| `--list` | Reading from a list file |
| `--table` (range) | CSV and TSV with CHROM/START/END columns |
| `--table` (position) | POS column with `--flank`, error without `--flank` |
| `--table` (extras) | Extra columns ignored, case-insensitive headers, any column order |
| `--sv-table` (range) | TSV with full range columns |
| `--sv-table` (position) | CSV with POS columns + `--flank`, error without `--flank` |
| `-o` / `--output` | Writing to a single output file |
| `--output-dir` | Per-region files, SV paired files |
| FAI handling | Auto-build, `--no-build-fai` error |
| Error cases | No regions, unknown contig, conflicting output flags, missing FASTA |

## Benchmarking

For performance benchmarks, use a real genome reference (e.g. GRCh38) and a
large region set. The `time` command gives a quick measurement:

```bash
# Build optimised binary
cargo build --release

# Benchmark with inline regions
time ./target/release/multiseqex GRCh38.fa \
    --regions chr1:1-1000000,chr2:1-1000000,chr3:1-1000000 \
    -o /dev/null

# Benchmark with a large table
time ./target/release/multiseqex GRCh38.fa \
    --table large_regions.tsv \
    -o /dev/null

# Compare single-threaded vs multi-threaded
time ./target/release/multiseqex GRCh38.fa --table large_regions.tsv --threads 1 -o /dev/null
time ./target/release/multiseqex GRCh38.fa --table large_regions.tsv --threads 4 -o /dev/null
time ./target/release/multiseqex GRCh38.fa --table large_regions.tsv -o /dev/null
```

### Comparison with samtools

To compare against `samtools faidx`:

```bash
# Extract the same regions with samtools
time samtools faidx GRCh38.fa chr1:1-1000000 chr2:1-1000000 chr3:1-1000000 > /dev/null

# For table input, convert to samtools format first
awk -F'\t' 'NR>1 {print $1":"$2"-"$3}' large_regions.tsv | \
    time xargs samtools faidx GRCh38.fa > /dev/null
```

`multiseqex` should outperform samtools when extracting many regions, thanks to
parallel I/O across multiple threads.

## Linting

```bash
cargo clippy
```

The project targets zero clippy warnings.
