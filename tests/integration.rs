use assert_cmd::cargo::cargo_bin_cmd;
use predicates::prelude::*;
use std::fs;
use tempfile::TempDir;

fn fixture(name: &str) -> String {
    format!("tests/fixtures/{name}")
}

fn cmd() -> assert_cmd::Command {
    cargo_bin_cmd!("multiseqex")
}

// ─── Basic extraction via --regions ──────────────────────────────────────────

#[test]
fn regions_single() {
    cmd()
        .arg(fixture("test.fa"))
        .args(["--regions", "chr1:1-10"])
        .assert()
        .success()
        .stdout(predicate::str::contains(">chr1:1-10"))
        .stdout(predicate::str::contains("AAACCCGGGT"));
}

#[test]
fn regions_multiple() {
    cmd()
        .arg(fixture("test.fa"))
        .args(["--regions", "chr1:1-10,chr2:1-10,chr3:1-10"])
        .assert()
        .success()
        .stdout(predicate::str::contains("AAACCCGGGT"))
        .stdout(predicate::str::contains("TTTTTTTTTT"))
        .stdout(predicate::str::contains("ATCGATCGAT"));
}

#[test]
fn regions_clamped_to_contig_length() {
    // chr3 is 72 bases. Requesting past the end should clamp.
    cmd()
        .arg(fixture("test.fa"))
        .args(["--regions", "chr3:60-999"])
        .assert()
        .success()
        .stdout(predicate::str::contains(">chr3:60-72"));
}

// ─── --list file ─────────────────────────────────────────────────────────────

#[test]
fn list_file() {
    cmd()
        .arg(fixture("test.fa"))
        .args(["--list", &fixture("regions.txt")])
        .assert()
        .success()
        .stdout(predicate::str::contains(">chr1:1-10"))
        .stdout(predicate::str::contains(">chr2:1-10"))
        .stdout(predicate::str::contains(">chr3:1-10"));
}

// ─── --table (CSV range mode) ────────────────────────────────────────────────

#[test]
fn table_csv_range() {
    cmd()
        .arg(fixture("test.fa"))
        .args(["--table", &fixture("table_range.csv")])
        .assert()
        .success()
        .stdout(predicate::str::contains(">chr1:1-10"))
        .stdout(predicate::str::contains("AAACCCGGGT"))
        .stdout(predicate::str::contains(">chr2:1-10"))
        .stdout(predicate::str::contains(">chr3:1-10"));
}

// ─── --table (TSV range mode with extra GENE column) ─────────────────────────

#[test]
fn table_tsv_range_extra_columns() {
    cmd()
        .arg(fixture("test.fa"))
        .args(["--table", &fixture("table_range.tsv")])
        .assert()
        .success()
        .stdout(predicate::str::contains(">chr1:1-10"))
        .stdout(predicate::str::contains(">chr2:1-10"))
        .stdout(predicate::str::contains(">chr3:1-10"));
}

// ─── --table (CSV position mode with --flank) ────────────────────────────────

#[test]
fn table_csv_position_with_flank() {
    cmd()
        .arg(fixture("test.fa"))
        .args(["--table", &fixture("table_pos.csv"), "--flank", "5"])
        .assert()
        .success()
        .stdout(predicate::str::contains(">regionA chr1:5-15"))
        .stdout(predicate::str::contains(">regionB chr2:5-15"));
}

#[test]
fn table_csv_position_without_flank_errors() {
    cmd()
        .arg(fixture("test.fa"))
        .args(["--table", &fixture("table_pos.csv")])
        .assert()
        .failure()
        .stderr(predicate::str::contains("--flank is required"));
}

// ─── --sv-table (TSV range mode) ────────────────────────────────────────────

#[test]
fn sv_table_tsv_range() {
    cmd()
        .arg(fixture("test.fa"))
        .args(["--sv-table", &fixture("sv_table_range.tsv")])
        .assert()
        .success()
        .stdout(predicate::str::contains(">SV001 chr1:1-10"))
        .stdout(predicate::str::contains(">SV001 chr2:1-10"));
}

// ─── --sv-table (CSV position mode with --flank) ─────────────────────────────

#[test]
fn sv_table_csv_position_with_flank() {
    cmd()
        .arg(fixture("test.fa"))
        .args(["--sv-table", &fixture("sv_table_pos.csv"), "--flank", "5"])
        .assert()
        .success()
        // sv_table_pos.csv has no NAME column, so headers are plain
        .stdout(predicate::str::contains(">chr1:5-15"))
        .stdout(predicate::str::contains(">chr2:5-15"));
}

#[test]
fn sv_table_csv_position_without_flank_errors() {
    cmd()
        .arg(fixture("test.fa"))
        .args(["--sv-table", &fixture("sv_table_pos.csv")])
        .assert()
        .failure()
        .stderr(predicate::str::contains("--flank is required"));
}

// ─── --output (single file) ─────────────────────────────────────────────────

#[test]
fn output_to_single_file() {
    let tmp = TempDir::new().unwrap();
    let out = tmp.path().join("out.fa");

    cmd()
        .arg(fixture("test.fa"))
        .args(["--regions", "chr1:1-10,chr2:1-10"])
        .args(["--output", out.to_str().unwrap()])
        .assert()
        .success()
        .stdout(predicate::str::is_empty());

    let content = fs::read_to_string(&out).unwrap();
    assert!(content.contains(">chr1:1-10"));
    assert!(content.contains(">chr2:1-10"));
    assert!(content.contains("AAACCCGGGT"));
    assert!(content.contains("TTTTTTTTTT"));
}

// ─── --output-dir (per-region files) ─────────────────────────────────────────

#[test]
fn output_dir_per_region() {
    let tmp = TempDir::new().unwrap();
    let dir = tmp.path().join("seqs");

    cmd()
        .arg(fixture("test.fa"))
        .args(["--regions", "chr1:1-10,chr2:1-10"])
        .args(["--output-dir", dir.to_str().unwrap()])
        .assert()
        .success();

    assert!(dir.join("chr1_1_10.fa").exists());
    assert!(dir.join("chr2_1_10.fa").exists());

    let c1 = fs::read_to_string(dir.join("chr1_1_10.fa")).unwrap();
    assert!(c1.contains("AAACCCGGGT"));
}

// ─── --output-dir with --sv-table (paired files) ─────────────────────────────

#[test]
fn output_dir_sv_paired() {
    let tmp = TempDir::new().unwrap();
    let dir = tmp.path().join("sv_seqs");

    cmd()
        .arg(fixture("test.fa"))
        .args(["--sv-table", &fixture("sv_table_range.tsv")])
        .args(["--output-dir", dir.to_str().unwrap()])
        .assert()
        .success();

    let entries: Vec<_> = fs::read_dir(&dir)
        .unwrap()
        .filter_map(|e| e.ok())
        .collect();
    assert_eq!(entries.len(), 1, "SV pair should produce exactly 1 file");

    let content = fs::read_to_string(entries[0].path()).unwrap();
    assert!(content.contains(">SV001 chr1:1-10"));
    assert!(content.contains(">SV001 chr2:1-10"));
}

// ─── Error cases ─────────────────────────────────────────────────────────────

#[test]
fn no_regions_errors() {
    cmd()
        .arg(fixture("test.fa"))
        .assert()
        .failure()
        .stderr(predicate::str::contains("No regions provided"));
}

#[test]
fn unknown_contig_errors() {
    cmd()
        .arg(fixture("test.fa"))
        .args(["--regions", "chrZ:1-10"])
        .assert()
        .failure()
        .stderr(predicate::str::contains("Contigs not found"));
}

#[test]
fn output_and_output_dir_conflict() {
    let tmp = TempDir::new().unwrap();
    cmd()
        .arg(fixture("test.fa"))
        .args(["--regions", "chr1:1-10"])
        .args(["--output", "out.fa"])
        .args(["--output-dir", tmp.path().to_str().unwrap()])
        .assert()
        .failure()
        .stderr(predicate::str::contains("Cannot use both"));
}

#[test]
fn missing_fasta_errors() {
    cmd()
        .arg("nonexistent.fa")
        .args(["--regions", "chr1:1-10", "--no-build-fai"])
        .assert()
        .failure();
}

// ─── FAI auto-build ──────────────────────────────────────────────────────────

#[test]
fn auto_builds_fai() {
    let tmp = TempDir::new().unwrap();
    let fasta = tmp.path().join("auto.fa");
    fs::write(
        &fasta,
        ">chr1\nAAAAAAAAAAAAAAAAAAAAA\n>chr2\nCCCCCCCCCCCCCCCCCCCC\n",
    )
    .unwrap();

    cmd()
        .arg(fasta.to_str().unwrap())
        .args(["--regions", "chr1:1-5"])
        .assert()
        .success()
        .stdout(predicate::str::contains("AAAAA"));

    let fai = tmp.path().join("auto.fa.fai");
    assert!(fai.exists(), ".fai should have been auto-built");
}

#[test]
fn no_build_fai_flag_errors_when_missing() {
    let tmp = TempDir::new().unwrap();
    let fasta = tmp.path().join("nobuild.fa");
    fs::write(&fasta, ">chr1\nAAAAA\n").unwrap();

    cmd()
        .arg(fasta.to_str().unwrap())
        .args(["--regions", "chr1:1-5", "--no-build-fai"])
        .assert()
        .failure()
        .stderr(predicate::str::contains("Missing index"));
}

// ─── Case-insensitive headers ────────────────────────────────────────────────

#[test]
fn table_case_insensitive_headers() {
    let tmp = TempDir::new().unwrap();
    let csv_path = tmp.path().join("lower.csv");
    fs::write(&csv_path, "chrom,start,end\nchr1,1,10\n").unwrap();

    cmd()
        .arg(fixture("test.fa"))
        .args(["--table", csv_path.to_str().unwrap()])
        .assert()
        .success()
        .stdout(predicate::str::contains(">chr1:1-10"));
}

// ─── Column order independence ───────────────────────────────────────────────

#[test]
fn table_columns_in_any_order() {
    let tmp = TempDir::new().unwrap();
    let csv_path = tmp.path().join("reordered.csv");
    // END before START, NAME first
    fs::write(&csv_path, "NAME,END,CHROM,START\nfoo,10,chr1,1\n").unwrap();

    cmd()
        .arg(fixture("test.fa"))
        .args(["--table", csv_path.to_str().unwrap()])
        .assert()
        .success()
        .stdout(predicate::str::contains(">foo chr1:1-10"))
        .stdout(predicate::str::contains("AAACCCGGGT"));
}

// ─── NAME in output headers ─────────────────────────────────────────────────

#[test]
fn name_absent_gives_plain_header() {
    // table_range.csv has no NAME column
    cmd()
        .arg(fixture("test.fa"))
        .args(["--table", &fixture("table_range.csv")])
        .assert()
        .success()
        .stdout(predicate::str::contains(">chr1:1-10"));
}

#[test]
fn name_present_appears_in_header() {
    let tmp = TempDir::new().unwrap();
    let csv_path = tmp.path().join("named.csv");
    fs::write(&csv_path, "CHROM,START,END,NAME\nchr1,1,10,myregion\n").unwrap();

    cmd()
        .arg(fixture("test.fa"))
        .args(["--table", csv_path.to_str().unwrap()])
        .assert()
        .success()
        .stdout(predicate::str::contains(">myregion chr1:1-10"));
}

// ─── NAME in per-region filenames ────────────────────────────────────────────

#[test]
fn output_dir_uses_name_in_filename() {
    let tmp = TempDir::new().unwrap();
    let csv_path = tmp.path().join("named.csv");
    fs::write(&csv_path, "CHROM,START,END,NAME\nchr1,1,10,myregion\n").unwrap();
    let dir = tmp.path().join("seqs");

    cmd()
        .arg(fixture("test.fa"))
        .args(["--table", csv_path.to_str().unwrap()])
        .args(["--output-dir", dir.to_str().unwrap()])
        .assert()
        .success();

    assert!(dir.join("myregion_1_10.fa").exists());
    let content = fs::read_to_string(dir.join("myregion_1_10.fa")).unwrap();
    assert!(content.contains(">myregion chr1:1-10"));
}

// ─── Inline position+flank syntax ───────────────────────────────────────────

#[test]
fn regions_position_plus_flank() {
    cmd()
        .arg(fixture("test.fa"))
        .args(["--regions", "chr1:10+5"])
        .assert()
        .success()
        .stdout(predicate::str::contains(">chr1:5-15"));
}

// ─── Combining multiple input sources ────────────────────────────────────────

#[test]
fn combine_regions_and_list() {
    cmd()
        .arg(fixture("test.fa"))
        .args(["--regions", "chr1:1-10"])
        .args(["--list", &fixture("regions.txt")])
        .assert()
        .success()
        // regions.txt has chr1:1-10, chr2:1-10, chr3:1-10
        .stdout(predicate::str::contains(">chr1:1-10"))
        .stdout(predicate::str::contains(">chr2:1-10"))
        .stdout(predicate::str::contains(">chr3:1-10"));
}

// ─── FASTA line wrapping ─────────────────────────────────────────────────────

#[test]
fn output_wraps_at_60_characters() {
    // Extract a region longer than 60bp to verify wrapping
    cmd()
        .arg(fixture("test.fa"))
        .args(["--regions", "chr1:1-80"])
        .assert()
        .success()
        .stdout(predicate::str::contains(">chr1:1-80\n"))
        // First wrapped line should be exactly 60 chars
        .stdout(predicate::function(|output: &str| {
            let lines: Vec<&str> = output.lines().collect();
            // line 0 is header, line 1 is first seq line
            lines.len() >= 3 && lines[1].len() == 60 && lines[2].len() == 20
        }));
}

// ─── Position mode with small position clamps start to 1 ────────────────────

#[test]
fn table_position_clamps_start_to_one() {
    let tmp = TempDir::new().unwrap();
    let csv_path = tmp.path().join("edge.csv");
    // POS=3 with flank=10 would give start=-7, should clamp to 1
    fs::write(&csv_path, "CHROM,POS\nchr1,3\n").unwrap();

    cmd()
        .arg(fixture("test.fa"))
        .args(["--table", csv_path.to_str().unwrap(), "--flank", "10"])
        .assert()
        .success()
        .stdout(predicate::str::contains(">chr1:1-13"));
}

// ─── Swapped start/end in region string ──────────────────────────────────────

#[test]
fn regions_swapped_start_end() {
    // start > end should be normalized
    cmd()
        .arg(fixture("test.fa"))
        .args(["--regions", "chr1:10-1"])
        .assert()
        .success()
        .stdout(predicate::str::contains(">chr1:1-10"))
        .stdout(predicate::str::contains("AAACCCGGGT"));
}

// ─── SV paired file uses NAME ────────────────────────────────────────────────

#[test]
fn output_dir_sv_paired_uses_name() {
    let tmp = TempDir::new().unwrap();
    let dir = tmp.path().join("sv_named");

    cmd()
        .arg(fixture("test.fa"))
        .args(["--sv-table", &fixture("sv_table_range.tsv")])
        .args(["--output-dir", dir.to_str().unwrap()])
        .assert()
        .success();

    // sv_table_range.tsv has NAME=SV001
    let expected_file = dir.join("SV001_chr1_1_10_chr2_1_10.fa");
    assert!(expected_file.exists(), "SV file should use NAME in filename");
    let content = fs::read_to_string(&expected_file).unwrap();
    assert!(content.contains(">SV001 chr1:1-10"));
    assert!(content.contains(">SV001 chr2:1-10"));
}
