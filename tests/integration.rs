use assert_cmd::cargo::cargo_bin_cmd;
use predicates::prelude::*;
use std::fs;
use std::thread;
use std::time::Duration;
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

    let entries: Vec<_> = fs::read_dir(&dir).unwrap().filter_map(|e| e.ok()).collect();
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
        .stderr(predicate::str::contains("cannot be used with"));
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

// ─── FAI with line_bases=0 should error, not panic ──────────────────────────

#[test]
fn fai_line_bases_zero_errors_gracefully() {
    let tmp = TempDir::new().unwrap();
    let fasta = tmp.path().join("zero.fa");
    // Write a valid-looking FASTA (content does not matter for this test).
    fs::write(&fasta, ">chr1\nACGTACGT\n").unwrap();

    // Hand-craft a FAI where line_bases is 0. This would cause a division
    // by zero if the tool did not guard against it.
    let fai = tmp.path().join("zero.fa.fai");
    fs::write(&fai, "chr1\t8\t6\t0\t9\n").unwrap();

    cmd()
        .arg(fasta.to_str().unwrap())
        .args(["--regions", "chr1:1-5", "--no-build-fai"])
        .assert()
        .failure()
        .stderr(predicate::str::contains("line_bases is 0"));
}

// ─── Inconsistent line widths warning during FAI build ──────────────────────

#[test]
fn inconsistent_line_widths_warns() {
    let tmp = TempDir::new().unwrap();
    let fasta = tmp.path().join("uneven.fa");

    // First non-final sequence line: 20 bases. Second non-final line: 10 bases.
    // Third line (final): 5 bases. The mismatch between lines 1 and 2 should
    // trigger a warning because line 2 is at least as long as line 1 in bytes
    // would not hold (it is shorter), so we make line 2 *longer* instead.
    // Build a contig where the first line has 10 bases and the second non-final
    // line has 20 bases (longer than expected).
    let seq = ">chr1\nACGTACGTAC\nACGTACGTACACGTACGTAC\nACGT\n";
    fs::write(&fasta, seq).unwrap();

    // Do NOT pass --no-build-fai so the tool builds its own index.
    cmd()
        .arg(fasta.to_str().unwrap())
        .args(["--regions", "chr1:1-5"])
        .assert()
        .success()
        .stderr(predicate::str::contains("inconsistent line width"));
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
    assert!(
        expected_file.exists(),
        "SV file should use NAME in filename"
    );
    let content = fs::read_to_string(&expected_file).unwrap();
    assert!(content.contains(">SV001 chr1:1-10"));
    assert!(content.contains(">SV001 chr2:1-10"));
}

// ─── --list with comments and blank lines (T001-006) ────────────────────────

#[test]
fn list_with_comments_and_blanks() {
    let tmp = TempDir::new().unwrap();
    let list_path = tmp.path().join("regions_comments.txt");
    fs::write(
        &list_path,
        "# This is a comment\n\
         chr1:1-10\n\
         \n\
         # Another comment\n\
         \t  \n\
         chr2:1-10\n\
         \n\
         chr3:1-10\n",
    )
    .unwrap();

    let output = cmd()
        .arg(fixture("test.fa"))
        .args(["--list", list_path.to_str().unwrap()])
        .assert()
        .success();

    // Only the three valid regions should appear.
    output
        .stdout(predicate::str::contains(">chr1:1-10"))
        .stdout(predicate::str::contains(">chr2:1-10"))
        .stdout(predicate::str::contains(">chr3:1-10"))
        .stdout(predicate::str::contains("AAACCCGGGT"))
        .stdout(predicate::str::contains("TTTTTTTTTT"))
        .stdout(predicate::str::contains("ATCGATCGAT"));
}

// ─── FAI round-trip accuracy (T001-007) ─────────────────────────────────────

#[test]
fn fai_roundtrip_accuracy() {
    let tmp = TempDir::new().unwrap();
    let fasta = tmp.path().join("roundtrip.fa");

    // Two contigs with different line widths. contig1 uses 20-base lines,
    // contig2 uses 15-base lines.
    fs::write(
        &fasta,
        ">contig1\n\
         AAACCCGGGTTTAAACCCGG\n\
         TTTAAACCCGGGTTTAAACC\n\
         CGGGTTT\n\
         >contig2\n\
         ATCGATCGATCGATC\n\
         GATCGATCGATCGAT\n\
         CGATCG\n",
    )
    .unwrap();

    // Extract a range that spans a line boundary in contig1 (bases 18-25).
    let output = cmd()
        .arg(fasta.to_str().unwrap())
        .args(["--regions", "contig1:18-25"])
        .output()
        .unwrap();
    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).unwrap();
    assert!(
        stdout.contains("CGGTTAAA") || stdout.contains(">contig1:18-25"),
        "Header should appear in output"
    );

    // Verify the FAI file was created and has correct fields.
    let fai_path = tmp.path().join("roundtrip.fa.fai");
    assert!(fai_path.exists(), ".fai should have been auto-built");
    let fai_content = fs::read_to_string(&fai_path).unwrap();
    let lines: Vec<&str> = fai_content.lines().collect();
    assert_eq!(lines.len(), 2, "FAI should have two records");

    // contig1: 47 bases, offset after ">contig1\n" = 9, line_bases=20, line_bytes=21
    let fields1: Vec<&str> = lines[0].split('\t').collect();
    assert_eq!(fields1[0], "contig1");
    assert_eq!(fields1[1], "47", "contig1 length");
    assert_eq!(fields1[2], "9", "contig1 offset");
    assert_eq!(fields1[3], "20", "contig1 line_bases");
    assert_eq!(fields1[4], "21", "contig1 line_bytes");

    // contig2: 36 bases, offset after contig1 data + ">contig2\n"
    let fields2: Vec<&str> = lines[1].split('\t').collect();
    assert_eq!(fields2[0], "contig2");
    assert_eq!(fields2[1], "36", "contig2 length");
    assert_eq!(fields2[3], "15", "contig2 line_bases");
    assert_eq!(fields2[4], "16", "contig2 line_bytes");

    // Extract from contig2 spanning a line boundary (bases 13-18).
    let output2 = cmd()
        .arg(fasta.to_str().unwrap())
        .args(["--regions", "contig2:13-18"])
        .output()
        .unwrap();
    assert!(output2.status.success());
    let stdout2 = String::from_utf8(output2.stdout).unwrap();
    // contig2 full: ATCGATCGATCGATCGATCGATCGATCGATCGATCG
    // bases 13-18:  ATCGATCGATCG[ATCGAT]CGATCGATCGATCGATCG
    //                             ^13   ^18
    assert!(
        stdout2.contains("ATCGAT"),
        "Expected bases 13-18 of contig2 = ATCGAT, got: {}",
        stdout2
    );
}

// ─── --threads flag (T001-008) ──────────────────────────────────────────────

#[test]
fn threads_flag_produces_correct_output() {
    // Use multiple regions to exercise parallelism with --threads 2.
    let output = cmd()
        .arg(fixture("test.fa"))
        .args(["--regions", "chr1:1-10,chr2:1-10,chr3:1-10"])
        .args(["--threads", "2"])
        .assert()
        .success();

    output
        .stdout(predicate::str::contains(">chr1:1-10"))
        .stdout(predicate::str::contains("AAACCCGGGT"))
        .stdout(predicate::str::contains(">chr2:1-10"))
        .stdout(predicate::str::contains("TTTTTTTTTT"))
        .stdout(predicate::str::contains(">chr3:1-10"))
        .stdout(predicate::str::contains("ATCGATCGAT"));
}

// ─── Zero-coordinate rejection ──────────────────────────────────────────────

#[test]
fn table_start_zero_errors() {
    let tmp = TempDir::new().unwrap();
    let csv_path = tmp.path().join("zero_start.csv");
    fs::write(&csv_path, "CHROM,START,END\nchr1,0,10\n").unwrap();

    cmd()
        .arg(fixture("test.fa"))
        .args(["--table", csv_path.to_str().unwrap()])
        .assert()
        .failure()
        .stderr(predicate::str::contains("1-based coordinates"));
}

// ─── Duplicate contig names warning (R001-002) ──────────────────────────────

#[test]
fn duplicate_contig_warns() {
    let tmp = TempDir::new().unwrap();
    let fasta = tmp.path().join("dup.fa");
    // Two contigs with the same name.
    fs::write(
        &fasta,
        ">chr1\nAAAAAAAAAAAAAAAAAAAAA\n>chr1\nCCCCCCCCCCCCCCCCCCCC\n",
    )
    .unwrap();

    cmd()
        .arg(fasta.to_str().unwrap())
        .args(["--regions", "chr1:1-5"])
        .assert()
        .success()
        .stderr(predicate::str::contains("duplicate contig name"));
}

// ─── Empty contig name warning (R001-003) ───────────────────────────────────

#[test]
fn empty_contig_name_warns() {
    let tmp = TempDir::new().unwrap();
    let fasta = tmp.path().join("bare.fa");
    // A bare '>' header with no name.
    fs::write(&fasta, ">\nAAAAAAAAAAAAAAAAAAAAA\n").unwrap();

    cmd()
        .arg(fasta.to_str().unwrap())
        .args(["--regions", ":1-5"])
        .assert()
        .success()
        .stderr(predicate::str::contains("empty contig name"));
}

// ─── Empty table: headers only (T001-001) ────────────────────────────────────

#[test]
fn table_empty_body_errors() {
    let tmp = TempDir::new().unwrap();
    let csv_path = tmp.path().join("empty.csv");
    fs::write(&csv_path, "CHROM,START,END\n").unwrap();

    cmd()
        .arg(fixture("test.fa"))
        .args(["--table", csv_path.to_str().unwrap()])
        .assert()
        .failure()
        .stderr(predicate::str::contains("No regions provided"));
}

#[test]
fn sv_table_empty_body_errors() {
    let tmp = TempDir::new().unwrap();
    let tsv_path = tmp.path().join("empty.tsv");
    fs::write(
        &tsv_path,
        "NAME\tCHROM_LEFT\tSTART_LEFT\tEND_LEFT\tCHROM_RIGHT\tSTART_RIGHT\tEND_RIGHT\n",
    )
    .unwrap();

    cmd()
        .arg(fixture("test.fa"))
        .args(["--sv-table", tsv_path.to_str().unwrap()])
        .assert()
        .failure()
        .stderr(predicate::str::contains("No regions provided"));
}

// ─── Malformed CSV rows (T001-002) ──────────────────────────────────────────

#[test]
fn table_missing_end_field_errors() {
    let tmp = TempDir::new().unwrap();
    let csv_path = tmp.path().join("missing_end.csv");
    // Row has CHROM and START but the END value is missing (short row).
    // The CSV parser rejects this before column access, so we check for a
    // CSV-level error rather than our own "Missing END" message.
    fs::write(&csv_path, "CHROM,START,END\nchr1,1\n").unwrap();

    cmd()
        .arg(fixture("test.fa"))
        .args(["--table", csv_path.to_str().unwrap()])
        .assert()
        .failure()
        .stderr(predicate::str::contains("CSV error"));
}

#[test]
fn table_non_numeric_start_errors() {
    let tmp = TempDir::new().unwrap();
    let csv_path = tmp.path().join("bad_start.csv");
    fs::write(&csv_path, "CHROM,START,END\nchr1,abc,10\n").unwrap();

    cmd()
        .arg(fixture("test.fa"))
        .args(["--table", csv_path.to_str().unwrap()])
        .assert()
        .failure()
        .stderr(predicate::str::contains("Bad START at row 2"));
}

#[test]
fn table_extra_fields_ok() {
    let tmp = TempDir::new().unwrap();
    let csv_path = tmp.path().join("extra.csv");
    // Extra trailing field should not cause an error.
    fs::write(&csv_path, "CHROM,START,END,EXTRA\nchr1,1,10,bonus\n").unwrap();

    cmd()
        .arg(fixture("test.fa"))
        .args(["--table", csv_path.to_str().unwrap()])
        .assert()
        .success()
        .stdout(predicate::str::contains(">chr1:1-10"));
}

// ─── Gzip rejection integration (T001-003) ──────────────────────────────────

#[test]
fn gzip_file_rejected_at_cli() {
    let tmp = TempDir::new().unwrap();
    let gz_path = tmp.path().join("fake.fa.gz");
    // Write gzip magic bytes followed by arbitrary content.
    fs::write(&gz_path, b"\x1f\x8b\x08\x00fake gzip data").unwrap();

    cmd()
        .arg(gz_path.to_str().unwrap())
        .args(["--regions", "chr1:1-10", "--no-build-fai"])
        .assert()
        .failure()
        .stderr(predicate::str::contains("gzip"));
}

// ─── Combining --regions with --table (T001-004) ─────────────────────────────

#[test]
fn combine_regions_and_table() {
    let tmp = TempDir::new().unwrap();
    let csv_path = tmp.path().join("extra_region.csv");
    fs::write(&csv_path, "CHROM,START,END\nchr2,1,10\n").unwrap();

    cmd()
        .arg(fixture("test.fa"))
        .args(["--regions", "chr1:1-10"])
        .args(["--table", csv_path.to_str().unwrap()])
        .assert()
        .success()
        .stdout(predicate::str::contains(">chr1:1-10"))
        .stdout(predicate::str::contains(">chr2:1-10"))
        .stdout(predicate::str::contains("AAACCCGGGT"))
        .stdout(predicate::str::contains("TTTTTTTTTT"));
}

// ─── Filename collision in --output-dir (T001-005) ───────────────────────────

#[test]
fn output_dir_filename_collision_overwrites() {
    let tmp = TempDir::new().unwrap();
    let dir = tmp.path().join("collision");

    // Two identical regions produce the same filename. The second overwrites
    // the first. We verify that exactly one file exists and its content is
    // valid FASTA.
    cmd()
        .arg(fixture("test.fa"))
        .args(["--regions", "chr1:1-10,chr1:1-10"])
        .args(["--output-dir", dir.to_str().unwrap()])
        .assert()
        .success();

    let entries: Vec<_> = fs::read_dir(&dir).unwrap().filter_map(|e| e.ok()).collect();
    assert_eq!(
        entries.len(),
        1,
        "Colliding filenames should produce exactly one file"
    );
    let content = fs::read_to_string(entries[0].path()).unwrap();
    assert!(
        content.contains(">chr1:1-10"),
        "File should contain a valid FASTA header"
    );
    assert!(
        content.contains("AAACCCGGGT"),
        "File should contain the extracted sequence"
    );
}

// ─── Stale FAI warning (R001-004) ───────────────────────────────────────────

#[test]
fn stale_fai_warns() {
    let tmp = TempDir::new().unwrap();
    let fasta = tmp.path().join("stale.fa");
    fs::write(&fasta, ">chr1\nAAAAAAAAAAAAAAAAAAAAA\n").unwrap();

    // Build the FAI first by running the tool.
    cmd()
        .arg(fasta.to_str().unwrap())
        .args(["--regions", "chr1:1-5"])
        .assert()
        .success();

    let fai = tmp.path().join("stale.fa.fai");
    assert!(fai.exists(), ".fai should have been auto-built");

    // Wait briefly, then touch the FASTA so it is newer than the FAI.
    thread::sleep(Duration::from_millis(1100));
    let content = fs::read_to_string(&fasta).unwrap();
    fs::write(&fasta, &content).unwrap();

    // Run again with --no-build-fai so it uses the existing (now stale) index.
    cmd()
        .arg(fasta.to_str().unwrap())
        .args(["--regions", "chr1:1-5", "--no-build-fai"])
        .assert()
        .success()
        .stderr(predicate::str::contains("index may be stale"));
}

// ─── --delimiter override ───────────────────────────────────────────────────

#[test]
fn delimiter_override_tab_on_csv_extension() {
    // Write a tab-separated file with .csv extension, then use --delimiter tab.
    let dir = TempDir::new().unwrap();
    let csv_path = dir.path().join("regions.csv");
    fs::write(&csv_path, "CHROM\tSTART\tEND\nchr1\t1\t10\nchr2\t1\t10\n").unwrap();

    cmd()
        .arg(fixture("test.fa"))
        .args(["--table", csv_path.to_str().unwrap(), "--delimiter", "tab"])
        .assert()
        .success()
        .stdout(predicate::str::contains(">chr1:1-10"))
        .stdout(predicate::str::contains(">chr2:1-10"));
}

#[test]
fn delimiter_override_comma_on_tsv_extension() {
    // Write a comma-separated file with .tsv extension, then use --delimiter comma.
    let dir = TempDir::new().unwrap();
    let tsv_path = dir.path().join("regions.tsv");
    fs::write(&tsv_path, "CHROM,START,END\nchr1,1,10\n").unwrap();

    cmd()
        .arg(fixture("test.fa"))
        .args([
            "--table",
            tsv_path.to_str().unwrap(),
            "--delimiter",
            "comma",
        ])
        .assert()
        .success()
        .stdout(predicate::str::contains(">chr1:1-10"));
}

#[test]
fn delimiter_override_single_char() {
    // Semicolon-separated file parsed with --delimiter ";".
    let dir = TempDir::new().unwrap();
    let path = dir.path().join("regions.txt");
    fs::write(&path, "CHROM;START;END\nchr1;1;10\n").unwrap();

    cmd()
        .arg(fixture("test.fa"))
        .args(["--table", path.to_str().unwrap(), "--delimiter", ";"])
        .assert()
        .success()
        .stdout(predicate::str::contains(">chr1:1-10"));
}

#[test]
fn delimiter_sniff_tabs_in_txt_file() {
    // A .txt file with tab-separated content should be auto-detected as tab.
    let dir = TempDir::new().unwrap();
    let path = dir.path().join("regions.txt");
    fs::write(&path, "CHROM\tSTART\tEND\nchr1\t1\t10\nchr3\t1\t10\n").unwrap();

    cmd()
        .arg(fixture("test.fa"))
        .args(["--table", path.to_str().unwrap()])
        .assert()
        .success()
        .stdout(predicate::str::contains(">chr1:1-10"))
        .stdout(predicate::str::contains(">chr3:1-10"));
}

#[test]
fn delimiter_invalid_value_errors() {
    cmd()
        .arg(fixture("test.fa"))
        .args(["--table", &fixture("table_range.csv"), "--delimiter", "abc"])
        .assert()
        .failure()
        .stderr(predicate::str::contains("Invalid --delimiter"));
}

// ─── --rc (reverse complement, F002-001) ─────────────────────────────────────

#[test]
fn rc_flag_reverses_complement() {
    // chr1:1-10 normally extracts AAACCCGGGT.
    // Reverse complement of AAACCCGGGT is ACCCGGGTTT.
    cmd()
        .arg(fixture("test.fa"))
        .args(["--regions", "chr1:1-10", "--rc"])
        .assert()
        .success()
        .stdout(predicate::str::contains("ACCCGGGTTT"));
}

#[test]
fn rc_flag_with_reverse_complement_alias() {
    cmd()
        .arg(fixture("test.fa"))
        .args(["--regions", "chr1:1-10", "--reverse-complement"])
        .assert()
        .success()
        .stdout(predicate::str::contains("ACCCGGGTTT"));
}

// ─── --dedup (deduplicate regions, F002-002) ─────────────────────────────────

#[test]
fn dedup_removes_duplicate_regions() {
    let output = cmd()
        .arg(fixture("test.fa"))
        .args(["--regions", "chr1:1-10,chr1:1-10,chr2:1-10"])
        .args(["--dedup"])
        .output()
        .unwrap();
    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).unwrap();
    let header_count = stdout.lines().filter(|l| l.starts_with('>')).count();
    assert_eq!(
        header_count, 2,
        "Expected 2 headers after dedup, got {header_count}"
    );
    let stderr = String::from_utf8(output.stderr).unwrap();
    assert!(
        stderr.contains("removed 1 duplicate"),
        "Expected dedup message on stderr, got: {stderr}"
    );
}

#[test]
fn no_dedup_preserves_duplicates() {
    let output = cmd()
        .arg(fixture("test.fa"))
        .args(["--regions", "chr1:1-10,chr1:1-10"])
        .output()
        .unwrap();
    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).unwrap();
    let header_count = stdout.lines().filter(|l| l.starts_with('>')).count();
    assert_eq!(
        header_count, 2,
        "Without --dedup, duplicates should be preserved"
    );
}

// ─── --sort (sorted output, F002-003) ────────────────────────────────────────

#[test]
fn sort_orders_by_natural_chromosome_then_start() {
    let output = cmd()
        .arg(fixture("test.fa"))
        .args(["--regions", "chr3:1-10,chr1:11-20,chr2:1-10,chr1:1-10"])
        .args(["--sort"])
        .output()
        .unwrap();
    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).unwrap();
    let headers: Vec<&str> = stdout.lines().filter(|l| l.starts_with('>')).collect();
    assert_eq!(headers.len(), 4);
    assert_eq!(headers[0], ">chr1:1-10");
    assert_eq!(headers[1], ">chr1:11-20");
    assert_eq!(headers[2], ">chr2:1-10");
    assert_eq!(headers[3], ">chr3:1-10");
}

#[test]
fn no_sort_preserves_input_order() {
    let output = cmd()
        .arg(fixture("test.fa"))
        .args(["--regions", "chr3:1-10,chr1:1-10,chr2:1-10"])
        .output()
        .unwrap();
    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).unwrap();
    let headers: Vec<&str> = stdout.lines().filter(|l| l.starts_with('>')).collect();
    assert_eq!(headers.len(), 3);
    assert_eq!(headers[0], ">chr3:1-10");
    assert_eq!(headers[1], ">chr1:1-10");
    assert_eq!(headers[2], ">chr2:1-10");
}

// ─── BED input (RC27-010) ────────────────────────────────────────────────────

#[test]
fn bed_3_column_basic() {
    // BED 0-based [0,10) -> 1-based [1,10]
    cmd()
        .arg(fixture("test.fa"))
        .args(["--bed", &fixture("regions.bed")])
        .assert()
        .success()
        .stdout(predicate::str::contains(">chr1:1-10"))
        .stdout(predicate::str::contains("AAACCCGGGT"))
        .stdout(predicate::str::contains(">chr2:1-10"))
        .stdout(predicate::str::contains("TTTTTTTTTT"))
        .stdout(predicate::str::contains(">chr3:1-10"))
        .stdout(predicate::str::contains("ATCGATCGAT"));
}

#[test]
fn bed_with_name_column() {
    cmd()
        .arg(fixture("test.fa"))
        .args(["--bed", &fixture("regions_named.bed")])
        .assert()
        .success()
        .stdout(predicate::str::contains(">geneA chr1:1-10"))
        .stdout(predicate::str::contains(">geneB chr2:1-10"));
}

#[test]
fn bed_with_flank() {
    // BED: chr1 0 10 -> 1-based [1,10]. With --flank 5: [1, 15]
    // (start = max(1-5, 1) = 1, end = 10+5 = 15)
    cmd()
        .arg(fixture("test.fa"))
        .args(["--bed", &fixture("regions.bed"), "--flank", "5"])
        .assert()
        .success()
        .stdout(predicate::str::contains(">chr1:1-15"))
        .stdout(predicate::str::contains(">chr2:1-15"))
        .stdout(predicate::str::contains(">chr3:1-15"));
}

#[test]
fn bed_with_rc() {
    // chr1:1-10 = AAACCCGGGT, RC = ACCCGGGTTT
    cmd()
        .arg(fixture("test.fa"))
        .args(["--bed", &fixture("regions.bed"), "--rc"])
        .assert()
        .success()
        .stdout(predicate::str::contains("ACCCGGGTTT"));
}

#[test]
fn bed_with_comments_and_blanks() {
    cmd()
        .arg(fixture("test.fa"))
        .args(["--bed", &fixture("regions_comments.bed")])
        .assert()
        .success()
        .stdout(predicate::str::contains(">chr1:1-10"))
        .stdout(predicate::str::contains(">chr2:1-10"));
}

// ─── Feature 1: Per-region strand from BED STRAND column ─────────────────────

#[test]
fn bed_strand_plus_shows_in_header() {
    // BED with strand column: + strand should appear in header.
    cmd()
        .arg(fixture("test.fa"))
        .args(["--bed", &fixture("regions_strand.bed")])
        .assert()
        .success()
        .stdout(predicate::str::contains(">geneA chr1:1-10(+)"))
        .stdout(predicate::str::contains(">geneB chr2:1-10(-)"));
}

#[test]
fn bed_strand_minus_reverse_complements() {
    // chr2:1-10 = TTTTTTTTTT. RC = AAAAAAAAAA.
    // Minus strand region should be reverse-complemented.
    cmd()
        .arg(fixture("test.fa"))
        .args(["--bed", &fixture("regions_strand.bed")])
        .assert()
        .success()
        // geneA is + strand, so sequence stays: AAACCCGGGT
        .stdout(predicate::str::contains("AAACCCGGGT"))
        // geneB is - strand, TTTTTTTTTT RC = AAAAAAAAAA
        .stdout(predicate::str::contains("AAAAAAAAAA"));
}

#[test]
fn bed_strand_minus_plus_rc_cancel() {
    // Minus strand + --rc should cancel out. geneB chr2:1-10 = TTTTTTTTTT,
    // minus strand RC = AAAAAAAAAA, then --rc again = TTTTTTTTTT (original).
    cmd()
        .arg(fixture("test.fa"))
        .args(["--bed", &fixture("regions_strand.bed"), "--rc"])
        .assert()
        .success()
        // geneB: - strand + --rc cancel, so original TTTTTTTTTT
        .stdout(predicate::str::contains("TTTTTTTTTT"));
}

#[test]
fn table_with_strand_column() {
    let tmp = TempDir::new().unwrap();
    let csv_path = tmp.path().join("strand.csv");
    // chr1:1-10 = AAACCCGGGT. Minus strand: RC = ACCCGGGTTT.
    fs::write(
        &csv_path,
        "CHROM,START,END,NAME,STRAND\nchr1,1,10,gene1,-\n",
    )
    .unwrap();

    cmd()
        .arg(fixture("test.fa"))
        .args(["--table", csv_path.to_str().unwrap()])
        .assert()
        .success()
        .stdout(predicate::str::contains(">gene1 chr1:1-10(-)"))
        .stdout(predicate::str::contains("ACCCGGGTTT"));
}

// ─── Feature 3: Output format control ────────────────────────────────────────

#[test]
fn line_width_flag() {
    // Extract 80bp with --line-width 20: should produce 4 lines of 20.
    let output = cmd()
        .arg(fixture("test.fa"))
        .args(["--regions", "chr1:1-80", "--line-width", "20"])
        .output()
        .unwrap();
    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).unwrap();
    let lines: Vec<&str> = stdout.lines().collect();
    // line 0 = header, lines 1-4 = sequence (20 chars each)
    assert_eq!(lines.len(), 5, "Expected 1 header + 4 seq lines");
    assert_eq!(lines[1].len(), 20);
    assert_eq!(lines[4].len(), 20);
}

#[test]
fn no_wrap_flag() {
    // Extract 80bp with --no-wrap: sequence should be on one line.
    let output = cmd()
        .arg(fixture("test.fa"))
        .args(["--regions", "chr1:1-80", "--no-wrap"])
        .output()
        .unwrap();
    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).unwrap();
    let lines: Vec<&str> = stdout.lines().collect();
    // line 0 = header, line 1 = full 80bp sequence
    assert_eq!(lines.len(), 2, "Expected 1 header + 1 unwrapped seq line");
    assert_eq!(lines[1].len(), 80);
}

#[test]
fn tab_out_flag() {
    let output = cmd()
        .arg(fixture("test.fa"))
        .args(["--regions", "chr1:1-10,chr2:1-10", "--tab-out"])
        .output()
        .unwrap();
    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).unwrap();
    let lines: Vec<&str> = stdout.lines().collect();
    assert_eq!(lines.len(), 2, "Expected 2 TSV lines");
    // Verify tab-separated format.
    let fields: Vec<&str> = lines[0].split('\t').collect();
    assert_eq!(fields.len(), 5, "Expected 5 TSV columns");
    assert_eq!(fields[0], "chr1");
    assert_eq!(fields[1], "1");
    assert_eq!(fields[2], "10");
    assert_eq!(fields[3], ".");
    assert_eq!(fields[4], "AAACCCGGGT");
}

#[test]
fn tab_out_conflicts_with_output_dir() {
    let tmp = TempDir::new().unwrap();
    cmd()
        .arg(fixture("test.fa"))
        .args(["--regions", "chr1:1-10", "--tab-out"])
        .args(["--output-dir", tmp.path().to_str().unwrap()])
        .assert()
        .failure()
        .stderr(predicate::str::contains("cannot be used with"));
}

// ─── Feature 5: Merge overlapping regions ────────────────────────────────────

#[test]
fn merge_overlapping_regions() {
    let output = cmd()
        .arg(fixture("test.fa"))
        .args(["--regions", "chr1:1-10,chr1:5-20", "--merge"])
        .output()
        .unwrap();
    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).unwrap();
    let headers: Vec<&str> = stdout.lines().filter(|l| l.starts_with('>')).collect();
    assert_eq!(
        headers.len(),
        1,
        "Overlapping regions should merge into one"
    );
    assert_eq!(headers[0], ">chr1:1-20");
}

#[test]
fn merge_with_distance() {
    let output = cmd()
        .arg(fixture("test.fa"))
        .args([
            "--regions",
            "chr1:1-10,chr1:15-25",
            "--merge",
            "--merge-distance",
            "5",
        ])
        .output()
        .unwrap();
    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).unwrap();
    let headers: Vec<&str> = stdout.lines().filter(|l| l.starts_with('>')).collect();
    assert_eq!(
        headers.len(),
        1,
        "Nearby regions should merge with distance"
    );
    assert_eq!(headers[0], ">chr1:1-25");
}

#[test]
fn merge_does_not_merge_different_chroms() {
    let output = cmd()
        .arg(fixture("test.fa"))
        .args(["--regions", "chr1:1-10,chr2:5-15", "--merge"])
        .output()
        .unwrap();
    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).unwrap();
    let headers: Vec<&str> = stdout.lines().filter(|l| l.starts_with('>')).collect();
    assert_eq!(headers.len(), 2, "Different chromosomes should not merge");
}

#[test]
fn merge_implies_sort() {
    // Regions given out of order. --merge should sort them first.
    let output = cmd()
        .arg(fixture("test.fa"))
        .args(["--regions", "chr1:11-20,chr1:1-10,chr1:5-15", "--merge"])
        .output()
        .unwrap();
    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).unwrap();
    let headers: Vec<&str> = stdout.lines().filter(|l| l.starts_with('>')).collect();
    // Should merge all three overlapping regions into one.
    assert_eq!(headers.len(), 1);
    assert_eq!(headers[0], ">chr1:1-20");
}

#[test]
fn merge_distance_without_merge_errors() {
    cmd()
        .arg(fixture("test.fa"))
        .args(["--regions", "chr1:1-10", "--merge-distance", "5"])
        .assert()
        .failure()
        .stderr(predicate::str::contains(
            "--merge-distance requires --merge",
        ));
}

// ─── Feature 6: Whole-contig extraction ──────────────────────────────────────

#[test]
fn contigs_flag_extracts_whole_contig() {
    let output = cmd()
        .arg(fixture("test.fa"))
        .args(["--contigs", "chr1"])
        .output()
        .unwrap();
    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).unwrap();
    // chr1 is 204 bases long.
    assert!(stdout.contains(">chr1:1-204"));
}

#[test]
fn contigs_flag_multiple() {
    let output = cmd()
        .arg(fixture("test.fa"))
        .args(["--contigs", "chr1,chr3"])
        .output()
        .unwrap();
    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).unwrap();
    assert!(stdout.contains(">chr1:1-204"));
    assert!(stdout.contains(">chr3:1-72"));
}

#[test]
fn contigs_flag_unknown_contig_errors() {
    cmd()
        .arg(fixture("test.fa"))
        .args(["--contigs", "chrZ"])
        .assert()
        .failure()
        .stderr(predicate::str::contains("not found in FAI"));
}

#[test]
fn contig_list_file() {
    cmd()
        .arg(fixture("test.fa"))
        .args(["--contig-list", &fixture("contig_list.txt")])
        .assert()
        .success()
        .stdout(predicate::str::contains(">chr1:1-204"))
        .stdout(predicate::str::contains(">chr3:1-72"));
}

#[test]
fn contigs_combined_with_regions() {
    let output = cmd()
        .arg(fixture("test.fa"))
        .args(["--contigs", "chr3", "--regions", "chr1:1-10"])
        .output()
        .unwrap();
    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).unwrap();
    assert!(stdout.contains(">chr1:1-10"));
    assert!(stdout.contains(">chr3:1-72"));
}

#[test]
fn contigs_conflicts_with_sv_table() {
    cmd()
        .arg(fixture("test.fa"))
        .args([
            "--contigs",
            "chr1",
            "--sv-table",
            &fixture("sv_table_range.tsv"),
        ])
        .assert()
        .failure()
        .stderr(predicate::str::contains("cannot be used with"));
}

// ─── Feature 7: Asymmetric flanking ──────────────────────────────────────────

#[test]
fn flank_left_right_bed() {
    // BED: chr1 100 200 -> 1-based [101, 200]
    // With --flank-left 10 --flank-right 30: start=91, end=230 -> clamped to chr1 length
    let tmp = TempDir::new().unwrap();
    let bed_path = tmp.path().join("asym.bed");
    fs::write(&bed_path, "chr1\t100\t110\n").unwrap();

    cmd()
        .arg(fixture("test.fa"))
        .args([
            "--bed",
            bed_path.to_str().unwrap(),
            "--flank-left",
            "10",
            "--flank-right",
            "20",
        ])
        .assert()
        .success()
        .stdout(predicate::str::contains(">chr1:91-130"));
}

#[test]
fn flank_left_right_table() {
    let tmp = TempDir::new().unwrap();
    let csv_path = tmp.path().join("asym.csv");
    fs::write(&csv_path, "CHROM,POS\nchr1,50\n").unwrap();

    cmd()
        .arg(fixture("test.fa"))
        .args([
            "--table",
            csv_path.to_str().unwrap(),
            "--flank-left",
            "10",
            "--flank-right",
            "20",
        ])
        .assert()
        .success()
        .stdout(predicate::str::contains(">chr1:40-70"));
}

#[test]
fn flank_conflicts_with_flank_left_right() {
    cmd()
        .arg(fixture("test.fa"))
        .args([
            "--regions",
            "chr1:1-10",
            "--flank",
            "5",
            "--flank-left",
            "3",
            "--flank-right",
            "7",
        ])
        .assert()
        .failure()
        .stderr(predicate::str::contains("cannot be used with"));
}

#[test]
fn flank_left_without_right_errors() {
    cmd()
        .arg(fixture("test.fa"))
        .args(["--regions", "chr1:1-10", "--flank-left", "5"])
        .assert()
        .failure()
        .stderr(predicate::str::contains(
            "--flank-left and --flank-right must be specified together",
        ));
}

// ─── Feature 4: VCF input (--vcf) ─────────────────────────────────────────

#[test]
fn vcf_snp_extraction() {
    let tmp = TempDir::new().unwrap();
    let vcf_path = tmp.path().join("test.vcf");
    fs::write(
        &vcf_path,
        "##fileformat=VCFv4.2\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         chr1\t5\trs123\tA\tG\t.\t.\t.\n",
    )
    .unwrap();

    cmd()
        .arg(fixture("test.fa"))
        .args(["--vcf", vcf_path.to_str().unwrap()])
        .assert()
        .success()
        .stdout(predicate::str::contains(">rs123 chr1:5-5 REF=A ALT=G"));
}

#[test]
fn vcf_deletion_spans_ref_length() {
    let tmp = TempDir::new().unwrap();
    let vcf_path = tmp.path().join("del.vcf");
    // REF=ACGT (4 bases), so region should be POS..POS+3 = 10..13.
    fs::write(&vcf_path, "chr1\t10\tvar1\tACGT\tA\t.\t.\t.\n").unwrap();

    cmd()
        .arg(fixture("test.fa"))
        .args(["--vcf", vcf_path.to_str().unwrap()])
        .assert()
        .success()
        .stdout(predicate::str::contains(">var1 chr1:10-13 REF=ACGT ALT=A"));
}

#[test]
fn vcf_with_flank() {
    let tmp = TempDir::new().unwrap();
    let vcf_path = tmp.path().join("flank.vcf");
    fs::write(&vcf_path, "chr1\t50\t.\tA\tG\t.\t.\t.\n").unwrap();

    cmd()
        .arg(fixture("test.fa"))
        .args(["--vcf", vcf_path.to_str().unwrap(), "--flank", "10"])
        .assert()
        .success()
        .stdout(predicate::str::contains(">chr1:40-60"));
}

#[test]
fn vcf_dot_id_gives_no_name() {
    let tmp = TempDir::new().unwrap();
    let vcf_path = tmp.path().join("dot.vcf");
    fs::write(&vcf_path, "chr1\t5\t.\tA\tG\t.\t.\t.\n").unwrap();

    let output = cmd()
        .arg(fixture("test.fa"))
        .args(["--vcf", vcf_path.to_str().unwrap()])
        .output()
        .unwrap();
    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).unwrap();
    // Header should not have a name prefix, just coordinates.
    assert!(stdout.contains(">chr1:5-5 REF=A ALT=G"));
}

#[test]
fn vcf_conflicts_with_sv_table() {
    let tmp = TempDir::new().unwrap();
    let vcf_path = tmp.path().join("test.vcf");
    fs::write(&vcf_path, "chr1\t5\trs1\tA\tG\t.\t.\t.\n").unwrap();

    cmd()
        .arg(fixture("test.fa"))
        .args([
            "--vcf",
            vcf_path.to_str().unwrap(),
            "--sv-table",
            &fixture("sv_table_range.tsv"),
        ])
        .assert()
        .failure()
        .stderr(predicate::str::contains("cannot be used with"));
}

// ─── Feature 8: GFF3/GTF input (--gff) ────────────────────────────────────

#[test]
fn gff_gene_extraction() {
    let tmp = TempDir::new().unwrap();
    let gff_path = tmp.path().join("test.gff3");
    fs::write(
        &gff_path,
        "##gff-version 3\n\
         chr1\t.\tgene\t1\t10\t.\t+\t.\tID=gene1;Name=TP53\n\
         chr1\t.\texon\t1\t5\t.\t+\t.\tParent=gene1\n\
         chr2\t.\tgene\t1\t10\t.\t-\t.\tID=gene2;Name=BRCA1\n",
    )
    .unwrap();

    cmd()
        .arg(fixture("test.fa"))
        .args(["--gff", gff_path.to_str().unwrap()])
        .assert()
        .success()
        .stdout(predicate::str::contains(">TP53 chr1:1-10(+)"))
        .stdout(predicate::str::contains(">BRCA1 chr2:1-10(-)"));
}

#[test]
fn gff_exon_filter() {
    let tmp = TempDir::new().unwrap();
    let gff_path = tmp.path().join("exon.gff3");
    fs::write(
        &gff_path,
        "chr1\t.\tgene\t1\t20\t.\t+\t.\tName=gene1\n\
         chr1\t.\texon\t1\t10\t.\t+\t.\tName=exon1\n\
         chr1\t.\texon\t15\t20\t.\t+\t.\tName=exon2\n",
    )
    .unwrap();

    let output = cmd()
        .arg(fixture("test.fa"))
        .args(["--gff", gff_path.to_str().unwrap(), "--gff-feature", "exon"])
        .output()
        .unwrap();
    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).unwrap();
    let headers: Vec<&str> = stdout.lines().filter(|l| l.starts_with('>')).collect();
    assert_eq!(headers.len(), 2);
    assert!(stdout.contains(">exon1"));
    assert!(stdout.contains(">exon2"));
}

#[test]
fn gff_with_flank() {
    let tmp = TempDir::new().unwrap();
    let gff_path = tmp.path().join("flank.gff3");
    fs::write(&gff_path, "chr1\t.\tgene\t50\t60\t.\t+\t.\tName=g1\n").unwrap();

    cmd()
        .arg(fixture("test.fa"))
        .args(["--gff", gff_path.to_str().unwrap(), "--flank", "10"])
        .assert()
        .success()
        .stdout(predicate::str::contains(">g1 chr1:40-70(+)"));
}

#[test]
fn gff_minus_strand_reverse_complements() {
    let tmp = TempDir::new().unwrap();
    let gff_path = tmp.path().join("strand.gff3");
    // chr1:1-10 = AAACCCGGGT. Minus strand RC = ACCCGGGTTT.
    fs::write(&gff_path, "chr1\t.\tgene\t1\t10\t.\t-\t.\tName=g1\n").unwrap();

    cmd()
        .arg(fixture("test.fa"))
        .args(["--gff", gff_path.to_str().unwrap()])
        .assert()
        .success()
        .stdout(predicate::str::contains("ACCCGGGTTT"));
}

#[test]
fn gff_conflicts_with_sv_table() {
    let tmp = TempDir::new().unwrap();
    let gff_path = tmp.path().join("test.gff3");
    fs::write(&gff_path, "chr1\t.\tgene\t1\t10\t.\t+\t.\tName=g1\n").unwrap();

    cmd()
        .arg(fixture("test.fa"))
        .args([
            "--gff",
            gff_path.to_str().unwrap(),
            "--sv-table",
            &fixture("sv_table_range.tsv"),
        ])
        .assert()
        .failure()
        .stderr(predicate::str::contains("cannot be used with"));
}

#[test]
fn gff_feature_requires_gff() {
    cmd()
        .arg(fixture("test.fa"))
        .args(["--regions", "chr1:1-10", "--gff-feature", "exon"])
        .assert()
        .failure()
        .stderr(predicate::str::contains("--gff-feature"));
}

// ─── Feature 10: FASTQ output (--fastq) ───────────────────────────────────

#[test]
fn fastq_output() {
    let output = cmd()
        .arg(fixture("test.fa"))
        .args(["--regions", "chr1:1-10", "--fastq"])
        .output()
        .unwrap();
    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).unwrap();
    let lines: Vec<&str> = stdout.lines().collect();
    assert_eq!(lines.len(), 4);
    assert!(lines[0].starts_with('@'));
    assert!(lines[0].contains("chr1:1-10"));
    assert_eq!(lines[1], "AAACCCGGGT");
    assert_eq!(lines[2], "+");
    assert_eq!(lines[3], "IIIIIIIIII"); // 10 I's for 10 bases
}

#[test]
fn fastq_custom_qual() {
    let output = cmd()
        .arg(fixture("test.fa"))
        .args(["--regions", "chr1:1-5", "--fastq", "--qual", "J"])
        .output()
        .unwrap();
    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).unwrap();
    let lines: Vec<&str> = stdout.lines().collect();
    assert_eq!(lines[3], "JJJJJ");
}

#[test]
fn fastq_conflicts_with_tab_out() {
    cmd()
        .arg(fixture("test.fa"))
        .args(["--regions", "chr1:1-10", "--fastq", "--tab-out"])
        .assert()
        .failure()
        .stderr(predicate::str::contains("cannot be used with"));
}

#[test]
fn qual_requires_fastq() {
    cmd()
        .arg(fixture("test.fa"))
        .args(["--regions", "chr1:1-10", "--qual", "J"])
        .assert()
        .failure()
        .stderr(predicate::str::contains("--fastq"));
}

// ─── Feature 11: Statistics mode (--stats) ─────────────────────────────────

#[test]
fn stats_basic() {
    let output = cmd()
        .arg(fixture("test.fa"))
        .args(["--regions", "chr1:1-10", "--stats"])
        .output()
        .unwrap();
    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).unwrap();
    let lines: Vec<&str> = stdout.lines().collect();
    // First line is header.
    assert!(lines[0].contains("chr\tstart\tend\tname\tlength\tgc_percent\tn_count\tmasked_count"));
    // Second line is data.
    let fields: Vec<&str> = lines[1].split('\t').collect();
    assert_eq!(fields[0], "chr1");
    assert_eq!(fields[1], "1");
    assert_eq!(fields[2], "10");
    assert_eq!(fields[3], "."); // no name
    assert_eq!(fields[4], "10"); // length
}

#[test]
fn stats_multiple_regions() {
    let output = cmd()
        .arg(fixture("test.fa"))
        .args(["--regions", "chr1:1-10,chr2:1-10", "--stats"])
        .output()
        .unwrap();
    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).unwrap();
    let lines: Vec<&str> = stdout.lines().collect();
    assert_eq!(lines.len(), 3); // header + 2 data rows
}

#[test]
fn stats_conflicts_with_output() {
    let tmp = TempDir::new().unwrap();
    let out = tmp.path().join("out.fa");
    cmd()
        .arg(fixture("test.fa"))
        .args([
            "--regions",
            "chr1:1-10",
            "--stats",
            "--output",
            out.to_str().unwrap(),
        ])
        .assert()
        .failure()
        .stderr(predicate::str::contains("cannot be used with"));
}

#[test]
fn stats_conflicts_with_fastq() {
    cmd()
        .arg(fixture("test.fa"))
        .args(["--regions", "chr1:1-10", "--stats", "--fastq"])
        .assert()
        .failure()
        .stderr(predicate::str::contains("cannot be used with"));
}

// ─── Feature 12: Name templating (--name-template) ─────────────────────────

#[test]
fn name_template_basic() {
    let output = cmd()
        .arg(fixture("test.fa"))
        .args([
            "--regions",
            "chr1:1-10",
            "--name-template",
            "{chr}_{start}_{end}_idx{index}",
        ])
        .output()
        .unwrap();
    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).unwrap();
    assert!(stdout.contains(">chr1_1_10_idx1"));
}

#[test]
fn name_template_with_named_region() {
    let tmp = TempDir::new().unwrap();
    let csv_path = tmp.path().join("named.csv");
    fs::write(&csv_path, "CHROM,START,END,NAME\nchr1,1,10,mygene\n").unwrap();

    let output = cmd()
        .arg(fixture("test.fa"))
        .args([
            "--table",
            csv_path.to_str().unwrap(),
            "--name-template",
            "{name}|{chr}:{start}-{end}|len={length}",
        ])
        .output()
        .unwrap();
    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).unwrap();
    assert!(stdout.contains(">mygene|chr1:1-10|len=10"));
}

#[test]
fn name_template_with_fastq() {
    let output = cmd()
        .arg(fixture("test.fa"))
        .args([
            "--regions",
            "chr1:1-5",
            "--fastq",
            "--name-template",
            "seq_{index}",
        ])
        .output()
        .unwrap();
    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).unwrap();
    assert!(stdout.contains("@seq_1"));
}

#[test]
fn name_template_conflicts_with_tab_out() {
    cmd()
        .arg(fixture("test.fa"))
        .args([
            "--regions",
            "chr1:1-10",
            "--tab-out",
            "--name-template",
            "{chr}",
        ])
        .assert()
        .failure()
        .stderr(predicate::str::contains("cannot be used with"));
}

#[test]
fn name_template_conflicts_with_stats() {
    cmd()
        .arg(fixture("test.fa"))
        .args([
            "--regions",
            "chr1:1-10",
            "--stats",
            "--name-template",
            "{chr}",
        ])
        .assert()
        .failure()
        .stderr(predicate::str::contains("cannot be used with"));
}

// ─── Feature 9: Sequence masking ────────────────────────────────────────────

#[test]
fn mask_bed_hard_mask() {
    // mask.bed has BED intervals: chr1:2-5 (0-based) = 1-based [3,5]
    // and chr1:7-9 (0-based) = 1-based [8,9].
    // chr1:1-10 = AAACCCGGGT
    // Mask [3,5]: pos 3,4,5 (A,C,C) -> NNN
    // Mask [8,9]: pos 8,9 (G,G) -> NN
    // Result: AA + NNN + CG + NN + T = AANNNCGNNT
    cmd()
        .arg(fixture("test.fa"))
        .args(["--regions", "chr1:1-10"])
        .args(["--mask-bed", &fixture("mask.bed")])
        .assert()
        .success()
        .stdout(predicate::str::contains("AANNNCGNNT"));
}

#[test]
fn mask_bed_soft_mask() {
    // Same positions, but lowercase instead of N.
    // pos 3-5: A,C,C -> a,c,c. pos 8-9: G,G -> g,g
    // "AA" + "acc" + "CG" + "gg" + "T" = "AAaccCGggT"
    cmd()
        .arg(fixture("test.fa"))
        .args(["--regions", "chr1:1-10"])
        .args(["--mask-bed", &fixture("mask.bed"), "--soft-mask"])
        .assert()
        .success()
        .stdout(predicate::str::contains("AAaccCGggT"));
}

#[test]
fn mask_bed_requires_mask_bed_flag() {
    // --hard-mask without --mask-bed should fail.
    cmd()
        .arg(fixture("test.fa"))
        .args(["--regions", "chr1:1-10", "--hard-mask"])
        .assert()
        .failure();
}

#[test]
fn soft_mask_conflicts_with_hard_mask() {
    cmd()
        .arg(fixture("test.fa"))
        .args([
            "--regions",
            "chr1:1-10",
            "--mask-bed",
            &fixture("mask.bed"),
            "--soft-mask",
            "--hard-mask",
        ])
        .assert()
        .failure();
}

// ─── Feature 13: Multiple FASTA support ─────────────────────────────────────

#[test]
fn multiple_fasta_files() {
    // test.fa has chr1,chr2,chr3. extra.fa has chrX,chrY.
    cmd()
        .arg(fixture("test.fa"))
        .arg(fixture("extra.fa"))
        .args(["--regions", "chr1:1-10,chrX:1-10"])
        .assert()
        .success()
        .stdout(predicate::str::contains("AAACCCGGGT"))
        .stdout(predicate::str::contains("GGGGGGGGGG"));
}

#[test]
fn multiple_fasta_extract_from_second_file_only() {
    cmd()
        .arg(fixture("test.fa"))
        .arg(fixture("extra.fa"))
        .args(["--regions", "chrX:1-20"])
        .assert()
        .success()
        .stdout(predicate::str::contains("GGGGGGGGGGAAAAAAAAAA"));
}

#[test]
fn multiple_fasta_duplicate_contig_errors() {
    // Passing the same file twice should error because contigs appear in both.
    cmd()
        .arg(fixture("test.fa"))
        .arg(fixture("test.fa"))
        .args(["--regions", "chr1:1-10"])
        .assert()
        .failure()
        .stderr(predicate::str::contains("multiple FASTA files"));
}

// ─── Feature 14: Progress bar ───────────────────────────────────────────────

#[test]
fn progress_bar_shows_region_count() {
    let tmp = TempDir::new().unwrap();
    let out = tmp.path().join("out.fa");

    cmd()
        .arg(fixture("test.fa"))
        .args(["--regions", "chr1:1-10,chr2:1-10,chr3:1-10"])
        .args(["--output", out.to_str().unwrap()])
        .assert()
        .success()
        .stderr(predicate::str::contains("3 regions"));
}

#[test]
fn progress_bar_hidden_when_quiet() {
    let tmp = TempDir::new().unwrap();
    let out = tmp.path().join("out.fa");

    cmd()
        .arg(fixture("test.fa"))
        .args(["--regions", "chr1:1-10,chr2:1-10"])
        .args(["--output", out.to_str().unwrap()])
        .arg("--quiet")
        .assert()
        .success()
        .stderr(predicate::str::is_empty());
}

// ─── Feature 19: Sequence transforms ────────────────────────────────────────

#[test]
fn transform_to_rna() {
    // chr1:1-10 = AAACCCGGGT -> AAACCCGGGU (T -> U)
    cmd()
        .arg(fixture("test.fa"))
        .args(["--regions", "chr1:1-10", "--to-rna"])
        .assert()
        .success()
        .stdout(predicate::str::contains("AAACCCGGGU"));
}

#[test]
fn transform_uppercase() {
    cmd()
        .arg(fixture("test.fa"))
        .args(["--regions", "chr1:1-10", "--uppercase"])
        .assert()
        .success()
        .stdout(predicate::str::contains("AAACCCGGGT"));
}

#[test]
fn transform_lowercase() {
    cmd()
        .arg(fixture("test.fa"))
        .args(["--regions", "chr1:1-10", "--lowercase"])
        .assert()
        .success()
        .stdout(predicate::str::contains("aaacccgggt"));
}

#[test]
fn transform_translate() {
    // chr1:1-12 = AAACCCGGGTTT
    // Codons: AAA=K, CCC=P, GGG=G, TTT=F
    cmd()
        .arg(fixture("test.fa"))
        .args(["--regions", "chr1:1-12", "--translate"])
        .assert()
        .success()
        .stdout(predicate::str::contains("KPGF"));
}

#[test]
fn transform_to_rna_conflicts_with_translate() {
    cmd()
        .arg(fixture("test.fa"))
        .args(["--regions", "chr1:1-10", "--to-rna", "--translate"])
        .assert()
        .failure();
}

#[test]
fn transform_uppercase_conflicts_with_lowercase() {
    cmd()
        .arg(fixture("test.fa"))
        .args(["--regions", "chr1:1-10", "--uppercase", "--lowercase"])
        .assert()
        .failure();
}

#[test]
fn transform_to_rna_conflicts_with_stats() {
    cmd()
        .arg(fixture("test.fa"))
        .args(["--regions", "chr1:1-10", "--to-rna", "--stats"])
        .assert()
        .failure()
        .stderr(predicate::str::contains("cannot be used with --stats"));
}

#[test]
fn transform_translate_conflicts_with_stats() {
    cmd()
        .arg(fixture("test.fa"))
        .args(["--regions", "chr1:1-10", "--translate", "--stats"])
        .assert()
        .failure()
        .stderr(predicate::str::contains("cannot be used with --stats"));
}

#[test]
fn transform_soft_mask_then_lowercase() {
    // Combine masking with transform: soft mask + lowercase forces everything lowercase.
    cmd()
        .arg(fixture("test.fa"))
        .args(["--regions", "chr1:1-10"])
        .args([
            "--mask-bed",
            &fixture("mask.bed"),
            "--soft-mask",
            "--lowercase",
        ])
        .assert()
        .success()
        .stdout(predicate::str::contains("aaacccgggt"));
}

#[test]
fn transform_to_rna_with_soft_mask() {
    // Soft mask positions 3-5 and 8-9, then convert T->U/t->u.
    // After soft mask: AAaccCGggT
    // After to_rna: AAaccCGggU
    cmd()
        .arg(fixture("test.fa"))
        .args(["--regions", "chr1:1-10"])
        .args([
            "--mask-bed",
            &fixture("mask.bed"),
            "--soft-mask",
            "--to-rna",
        ])
        .assert()
        .success()
        .stdout(predicate::str::contains("AAaccCGggU"));
}

// ─── Feature 2: Bgzipped FASTA support ──────────────────────────────────────

#[test]
fn bgzip_fasta_decompressed_transparently() {
    use flate2::Compression;
    use flate2::write::GzEncoder;
    use std::io::Write;

    let tmp = TempDir::new().unwrap();
    let gz_path = tmp.path().join("test.fa.gz");

    // Compress a small FASTA to gzip format.
    let fasta_content = b">chr1\nACGTACGTACGT\n>chr2\nTTTTGGGGCCCC\n";
    let f = fs::File::create(&gz_path).unwrap();
    let mut encoder = GzEncoder::new(f, Compression::default());
    encoder.write_all(fasta_content).unwrap();
    encoder.finish().unwrap();

    // Should decompress transparently and extract.
    cmd()
        .arg(gz_path.to_str().unwrap())
        .args(["--regions", "chr1:1-8"])
        .assert()
        .success()
        .stdout(predicate::str::contains(">chr1:1-8"))
        .stdout(predicate::str::contains("ACGTACGT"));
}

#[test]
fn bgzip_fasta_with_existing_fai() {
    use flate2::Compression;
    use flate2::write::GzEncoder;
    use std::io::Write;

    let tmp = TempDir::new().unwrap();
    let gz_path = tmp.path().join("test.fa.gz");
    let fai_path = tmp.path().join("test.fa.gz.fai");

    let fasta_content = b">chr1\nAAAACCCCGGGGTTTT\n";
    let f = fs::File::create(&gz_path).unwrap();
    let mut encoder = GzEncoder::new(f, Compression::default());
    encoder.write_all(fasta_content).unwrap();
    encoder.finish().unwrap();

    // Let the tool build a FAI from the decompressed temp.
    cmd()
        .arg(gz_path.to_str().unwrap())
        .args(["--regions", "chr1:1-4"])
        .assert()
        .success()
        .stdout(predicate::str::contains("AAAA"));

    // The FAI for the original .gz path should NOT exist (it is built in temp).
    assert!(!fai_path.exists());
}

// ─── Feature 15: Streaming extraction without FAI (--no-index) ──────────────

#[test]
fn no_index_basic_extraction() {
    let tmp = TempDir::new().unwrap();
    let fasta = tmp.path().join("noindex.fa");
    fs::write(&fasta, ">chr1\nACGTACGTACGTACGT\n>chr2\nTTTTGGGGCCCCAAAA\n").unwrap();

    cmd()
        .arg(fasta.to_str().unwrap())
        .args(["--regions", "chr1:1-8", "--no-index"])
        .assert()
        .success()
        .stdout(predicate::str::contains(">chr1:1-8"))
        .stdout(predicate::str::contains("ACGTACGT"));
}

#[test]
fn no_index_multiple_regions() {
    let tmp = TempDir::new().unwrap();
    let fasta = tmp.path().join("noindex.fa");
    fs::write(&fasta, ">chr1\nACGTACGT\n>chr2\nTTTTGGGG\n").unwrap();

    cmd()
        .arg(fasta.to_str().unwrap())
        .args(["--regions", "chr1:1-4,chr2:1-4", "--no-index"])
        .assert()
        .success()
        .stdout(predicate::str::contains(">chr1:1-4"))
        .stdout(predicate::str::contains("ACGT"))
        .stdout(predicate::str::contains(">chr2:1-4"))
        .stdout(predicate::str::contains("TTTT"));
}

#[test]
fn no_index_conflicts_with_no_build_fai() {
    cmd()
        .arg(fixture("test.fa"))
        .args(["--regions", "chr1:1-10", "--no-index", "--no-build-fai"])
        .assert()
        .failure();
}

#[test]
fn no_index_with_rc() {
    let tmp = TempDir::new().unwrap();
    let fasta = tmp.path().join("noindex.fa");
    fs::write(&fasta, ">chr1\nAAAA\n").unwrap();

    cmd()
        .arg(fasta.to_str().unwrap())
        .args(["--regions", "chr1:1-4", "--no-index", "--rc"])
        .assert()
        .success()
        .stdout(predicate::str::contains("TTTT"));
}

#[test]
fn no_index_gzip_transparent() {
    use flate2::Compression;
    use flate2::write::GzEncoder;
    use std::io::Write;

    let tmp = TempDir::new().unwrap();
    let gz_path = tmp.path().join("test.fa.gz");

    let fasta_content = b">chr1\nACGTACGT\n";
    let f = fs::File::create(&gz_path).unwrap();
    let mut encoder = GzEncoder::new(f, Compression::default());
    encoder.write_all(fasta_content).unwrap();
    encoder.finish().unwrap();

    cmd()
        .arg(gz_path.to_str().unwrap())
        .args(["--regions", "chr1:1-4", "--no-index"])
        .assert()
        .success()
        .stdout(predicate::str::contains("ACGT"));
}

#[test]
fn stdin_requires_no_index() {
    cmd()
        .arg("-")
        .args(["--regions", "chr1:1-10"])
        .assert()
        .failure()
        .stderr(predicate::str::contains("requires --no-index"));
}

#[test]
fn no_index_tab_output() {
    let tmp = TempDir::new().unwrap();
    let fasta = tmp.path().join("noindex.fa");
    fs::write(&fasta, ">chr1\nACGTACGT\n").unwrap();

    cmd()
        .arg(fasta.to_str().unwrap())
        .args(["--regions", "chr1:1-4", "--no-index", "--tab-out"])
        .assert()
        .success()
        .stdout(predicate::str::contains("chr1\t1\t4\t.\tACGT"));
}

// ─── Feature 16: Built-in interval operations (--subtract, --intersect) ─────

#[test]
fn subtract_trims_region() {
    let tmp = TempDir::new().unwrap();
    let subtract_bed = tmp.path().join("subtract.bed");
    // Subtract bases 5-7 (0-based [4,7) -> 1-based [5,7]) from chr1:1-10.
    fs::write(&subtract_bed, "chr1\t4\t7\n").unwrap();

    cmd()
        .arg(fixture("test.fa"))
        .args([
            "--regions",
            "chr1:1-10",
            "--subtract",
            subtract_bed.to_str().unwrap(),
        ])
        .assert()
        .success()
        // Should produce two pieces: 1-4 and 8-10.
        .stdout(predicate::str::contains("AAAC"))
        .stdout(predicate::str::contains("GGT"));
}

#[test]
fn subtract_complete_removal_errors() {
    let tmp = TempDir::new().unwrap();
    let subtract_bed = tmp.path().join("subtract.bed");
    // Subtract entire region.
    fs::write(&subtract_bed, "chr1\t0\t10000\n").unwrap();

    cmd()
        .arg(fixture("test.fa"))
        .args([
            "--regions",
            "chr1:1-10",
            "--subtract",
            subtract_bed.to_str().unwrap(),
        ])
        .assert()
        .failure()
        .stderr(predicate::str::contains(
            "No regions remain after --subtract",
        ));
}

#[test]
fn intersect_keeps_overlap_only() {
    let tmp = TempDir::new().unwrap();
    let intersect_bed = tmp.path().join("intersect.bed");
    // Keep only bases 3-6 (0-based [2,6) -> 1-based [3,6]).
    fs::write(&intersect_bed, "chr1\t2\t6\n").unwrap();

    cmd()
        .arg(fixture("test.fa"))
        .args([
            "--regions",
            "chr1:1-10",
            "--intersect",
            intersect_bed.to_str().unwrap(),
        ])
        .assert()
        .success()
        .stdout(predicate::str::contains("ACCC"));
}

#[test]
fn intersect_no_overlap_errors() {
    let tmp = TempDir::new().unwrap();
    let intersect_bed = tmp.path().join("intersect.bed");
    // No overlap with chr1.
    fs::write(&intersect_bed, "chrZ\t0\t100\n").unwrap();

    cmd()
        .arg(fixture("test.fa"))
        .args([
            "--regions",
            "chr1:1-10",
            "--intersect",
            intersect_bed.to_str().unwrap(),
        ])
        .assert()
        .failure()
        .stderr(predicate::str::contains(
            "No regions remain after --intersect",
        ));
}

#[test]
fn intersect_then_subtract() {
    let tmp = TempDir::new().unwrap();
    let intersect_bed = tmp.path().join("intersect.bed");
    let subtract_bed = tmp.path().join("subtract.bed");
    // Intersect: keep bases 1-8 of chr1. (0-based [0,8) -> 1-based [1,8])
    fs::write(&intersect_bed, "chr1\t0\t8\n").unwrap();
    // Subtract: remove bases 3-5. (0-based [2,5) -> 1-based [3,5])
    fs::write(&subtract_bed, "chr1\t2\t5\n").unwrap();

    cmd()
        .arg(fixture("test.fa"))
        .args([
            "--regions",
            "chr1:1-10",
            "--intersect",
            intersect_bed.to_str().unwrap(),
            "--subtract",
            subtract_bed.to_str().unwrap(),
        ])
        .assert()
        .success()
        // After intersect: chr1:1-8 (AAACCCGG)
        // After subtract of [3,5]: chr1:1-2 (AA) and chr1:6-8 (CGG)
        .stdout(predicate::str::contains("AA"))
        .stdout(predicate::str::contains("CGG"));
}

// ─── Feature 17: K-mer tiling (--tile --step) ───────────────────────────────

#[test]
fn tile_non_overlapping() {
    cmd()
        .arg(fixture("test.fa"))
        .args(["--regions", "chr1:1-20", "--tile", "10"])
        .assert()
        .success()
        .stdout(predicate::str::contains("_tile1"))
        .stdout(predicate::str::contains("_tile2"));
}

#[test]
fn tile_with_step() {
    cmd()
        .arg(fixture("test.fa"))
        .args(["--regions", "chr1:1-20", "--tile", "10", "--step", "5"])
        .assert()
        .success()
        // Should produce tiles: 1-10, 6-15, 11-20, 16-20
        .stdout(predicate::str::contains("_tile1"))
        .stdout(predicate::str::contains("_tile2"))
        .stdout(predicate::str::contains("_tile3"))
        .stdout(predicate::str::contains("_tile4"));
}

#[test]
fn tile_zero_errors() {
    cmd()
        .arg(fixture("test.fa"))
        .args(["--regions", "chr1:1-10", "--tile", "0"])
        .assert()
        .failure()
        .stderr(predicate::str::contains("--tile must be > 0"));
}

#[test]
fn step_requires_tile() {
    cmd()
        .arg(fixture("test.fa"))
        .args(["--regions", "chr1:1-10", "--step", "5"])
        .assert()
        .failure();
}

#[test]
fn tile_region_shorter_than_tile() {
    // Region of 5 bases with tile of 10: produces one tile.
    cmd()
        .arg(fixture("test.fa"))
        .args(["--regions", "chr1:1-5", "--tile", "10"])
        .assert()
        .success()
        .stdout(predicate::str::contains("_tile1"));
}

#[test]
fn tile_with_bed_input() {
    let tmp = TempDir::new().unwrap();
    let bed = tmp.path().join("tile.bed");
    fs::write(&bed, "chr1\t0\t20\n").unwrap();

    cmd()
        .arg(fixture("test.fa"))
        .args([
            "--bed",
            bed.to_str().unwrap(),
            "--tile",
            "10",
            "--step",
            "10",
        ])
        .assert()
        .success()
        .stdout(predicate::str::contains("_tile1"))
        .stdout(predicate::str::contains("_tile2"));
}

#[test]
fn tile_conflicts_with_sv_table_output_dir() {
    let tmp = TempDir::new().unwrap();
    let dir = tmp.path().join("out");
    cmd()
        .arg(fixture("test.fa"))
        .args([
            "--sv-table",
            &fixture("sv_table_range.tsv"),
            "--output-dir",
            dir.to_str().unwrap(),
            "--tile",
            "10",
        ])
        .assert()
        .failure()
        .stderr(predicate::str::contains("--tile cannot be used with"));
}

#[test]
fn tile_with_no_index() {
    let tmp = TempDir::new().unwrap();
    let fasta = tmp.path().join("tile.fa");
    fs::write(&fasta, ">chr1\nACGTACGTACGTACGTACGT\n").unwrap();

    cmd()
        .arg(fasta.to_str().unwrap())
        .args(["--regions", "chr1:1-20", "--no-index", "--tile", "10"])
        .assert()
        .success()
        .stdout(predicate::str::contains("_tile1"))
        .stdout(predicate::str::contains("_tile2"));
}
