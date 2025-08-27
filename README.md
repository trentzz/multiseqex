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
      --table <TABLE>            CSV/TSV with 3 columns: chr,start,end (header allowed) Delimiter auto-detected from extension: .tsv => tab, else comma
      --flank <FLANK>            Flank size (only used if table has chr,pos with 2 columns) [default: 0]
  -o, --output <OUTPUT>          Output FASTA file (single combined; default: stdout)
      --output-dir <OUTPUT_DIR>  Output directory for per-sequence FASTA files
      --threads <THREADS>        Number of worker threads (default: Rayon default, usually #cpus)
      --no-build-fai             Do not build a .fai if missing (error instead)
  -h, --help                     Print help
  -V, --version                  Print version
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
