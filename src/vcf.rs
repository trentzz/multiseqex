//! Minimal VCF parser for extracting regions from VCF files.
//!
//! Parses standard VCF (v4.x) format. Lines beginning with `#` are skipped.
//! Data lines are tab-separated with fields: CHROM, POS, ID, REF, ALT, ...
//! VCF coordinates are 1-based. Each record produces a region spanning
//! `POS..POS+len(REF)-1`, optionally extended by flanking.

use anyhow::{Context, Result, anyhow};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use crate::region::{Region, resolve_flanks};

/// A parsed VCF record with the fields relevant to region extraction.
#[derive(Debug, Clone)]
pub struct VcfRecord {
    pub chrom: String,
    pub pos: u64,
    pub id: Option<String>,
    pub ref_allele: String,
    pub alt_allele: String,
}

/// Format a FASTA header description from VCF REF and ALT fields.
pub fn vcf_description(rec: &VcfRecord) -> String {
    format!("REF={} ALT={}", rec.ref_allele, rec.alt_allele)
}

/// Parse a VCF file and return regions with optional flanking.
///
/// Each record produces a region from `POS` to `POS + len(REF) - 1`.
/// If flanking is specified, it is applied around that region.
/// The ID field (column 3) is used as the region name unless it is `.`.
pub fn parse_regions_vcf(
    path: &Path,
    flank: Option<u64>,
    flank_left: Option<u64>,
    flank_right: Option<u64>,
) -> Result<Vec<(Region, VcfRecord)>> {
    let f =
        File::open(path).with_context(|| format!("Cannot open VCF file: {}", path.display()))?;
    let reader = BufReader::new(f);
    let (fl, fr) = resolve_flanks(flank, flank_left, flank_right);
    let has_flank = flank.is_some() || flank_left.is_some() || flank_right.is_some();
    let mut results = Vec::new();

    for (line_num, line_result) in reader.lines().enumerate() {
        let line = line_result
            .with_context(|| format!("I/O error reading VCF file: {}", path.display()))?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = trimmed.split('\t').collect();
        if fields.len() < 5 {
            return Err(anyhow!(
                "VCF line {} has fewer than 5 fields: {}",
                line_num + 1,
                trimmed
            ));
        }
        let chrom = fields[0].to_string();
        let pos: u64 = fields[1]
            .parse()
            .with_context(|| format!("Bad POS at VCF line {}: {}", line_num + 1, fields[1]))?;
        if pos == 0 {
            return Err(anyhow!(
                "VCF line {} has POS=0; VCF uses 1-based coordinates",
                line_num + 1
            ));
        }
        let id_field = fields[2];
        let id = if id_field == "." {
            None
        } else {
            Some(id_field.to_string())
        };
        let ref_allele = fields[3].to_string();
        let alt_allele = fields[4].to_string();

        let ref_len = ref_allele.len() as u64;
        let region_start = pos;
        let region_end = pos.saturating_add(ref_len.saturating_sub(1));

        let (final_start, final_end) = if has_flank {
            (
                region_start.saturating_sub(fl).max(1),
                region_end.saturating_add(fr),
            )
        } else {
            (region_start, region_end)
        };

        let rec = VcfRecord {
            chrom: chrom.clone(),
            pos,
            id: id.clone(),
            ref_allele,
            alt_allele,
        };

        results.push((
            Region {
                name: id,
                chr: chrom,
                start: final_start,
                end: final_end,
                strand: None,
                alt_info: None,
            },
            rec,
        ));
    }
    Ok(results)
}

/// Expand VCF results for `--alt-seq` mode.
///
/// For each VCF record, produces one `Region` per ALT allele (splitting
/// multi-allelic sites on commas). Each region carries an `AltInfo` so the
/// substitution can be applied after extraction.
///
/// When `both` is true, a reference-sequence region (without `AltInfo`) is
/// emitted before each alt region.
pub fn expand_vcf_alt_seq(
    vcf_results: Vec<(Region, VcfRecord)>,
    both: bool,
) -> Vec<(Region, VcfRecord)> {
    use crate::region::AltInfo;

    let mut out = Vec::new();
    for (region, rec) in vcf_results {
        let alts: Vec<&str> = rec.alt_allele.split(',').collect();
        for alt in &alts {
            let alt = alt.trim();
            if alt.is_empty() || alt == "." {
                continue;
            }

            // Build a VcfRecord for this single ALT allele.
            let single_rec = VcfRecord {
                chrom: rec.chrom.clone(),
                pos: rec.pos,
                id: rec.id.clone(),
                ref_allele: rec.ref_allele.clone(),
                alt_allele: alt.to_string(),
            };

            if both {
                // Emit the reference sequence entry first (no AltInfo).
                let mut ref_region = region.clone();
                ref_region.alt_info = None;
                out.push((ref_region, single_rec.clone()));
            }

            // Emit the alt-seq entry.
            let mut alt_region = region.clone();
            alt_region.alt_info = Some(AltInfo {
                ref_allele: rec.ref_allele.clone(),
                alt_allele: alt.to_string(),
                variant_pos: rec.pos,
            });
            out.push((alt_region, single_rec));
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_snp_no_flank() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        std::io::Write::write_all(
            &mut tmp.as_file(),
            b"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t100\trs123\tA\tG\t.\t.\t.\n",
        )
        .unwrap();
        let results = parse_regions_vcf(tmp.path(), None, None, None).unwrap();
        assert_eq!(results.len(), 1);
        let (region, rec) = &results[0];
        assert_eq!(region.chr, "chr1");
        assert_eq!(region.start, 100);
        assert_eq!(region.end, 100); // SNP: len(REF)=1, so end=100
        assert_eq!(region.name.as_deref(), Some("rs123"));
        assert_eq!(rec.ref_allele, "A");
        assert_eq!(rec.alt_allele, "G");
    }

    #[test]
    fn parse_deletion_no_flank() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        std::io::Write::write_all(&mut tmp.as_file(), b"chr1\t50\t.\tACGT\tA\t.\t.\t.\n").unwrap();
        let results = parse_regions_vcf(tmp.path(), None, None, None).unwrap();
        assert_eq!(results.len(), 1);
        let (region, _) = &results[0];
        assert_eq!(region.start, 50);
        assert_eq!(region.end, 53); // len("ACGT")=4, so 50+3=53
        assert!(region.name.is_none()); // ID is "."
    }

    #[test]
    fn parse_with_symmetric_flank() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        std::io::Write::write_all(&mut tmp.as_file(), b"chr1\t100\trs1\tA\tG\t.\t.\t.\n").unwrap();
        let results = parse_regions_vcf(tmp.path(), Some(10), None, None).unwrap();
        let (region, _) = &results[0];
        assert_eq!(region.start, 90); // 100-10
        assert_eq!(region.end, 110); // 100+10
    }

    #[test]
    fn parse_with_asymmetric_flank() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        std::io::Write::write_all(&mut tmp.as_file(), b"chr1\t100\trs1\tA\tG\t.\t.\t.\n").unwrap();
        let results = parse_regions_vcf(tmp.path(), None, Some(5), Some(20)).unwrap();
        let (region, _) = &results[0];
        assert_eq!(region.start, 95); // 100-5
        assert_eq!(region.end, 120); // 100+20
    }

    #[test]
    fn flank_clamps_start_to_one() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        std::io::Write::write_all(&mut tmp.as_file(), b"chr1\t3\trs1\tA\tG\t.\t.\t.\n").unwrap();
        let results = parse_regions_vcf(tmp.path(), Some(10), None, None).unwrap();
        let (region, _) = &results[0];
        assert_eq!(region.start, 1);
        assert_eq!(region.end, 13);
    }

    #[test]
    fn vcf_description_format() {
        let rec = VcfRecord {
            chrom: "chr1".into(),
            pos: 100,
            id: Some("rs123".into()),
            ref_allele: "A".into(),
            alt_allele: "G".into(),
        };
        assert_eq!(vcf_description(&rec), "REF=A ALT=G");
    }

    #[test]
    fn skip_comment_and_header_lines() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        std::io::Write::write_all(
            &mut tmp.as_file(),
            b"##fileformat=VCFv4.2\n##INFO=<ID=DP>\n#CHROM\tPOS\tID\tREF\tALT\n\nchr1\t10\t.\tA\tG\t.\t.\t.\n",
        )
        .unwrap();
        let results = parse_regions_vcf(tmp.path(), None, None, None).unwrap();
        assert_eq!(results.len(), 1);
    }

    #[test]
    fn pos_zero_errors() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        std::io::Write::write_all(&mut tmp.as_file(), b"chr1\t0\trs1\tA\tG\t.\t.\t.\n").unwrap();
        let err = parse_regions_vcf(tmp.path(), None, None, None).unwrap_err();
        assert!(err.to_string().contains("1-based"));
    }

    #[test]
    fn expand_alt_seq_single_allele() {
        let results = vec![(
            Region {
                name: Some("rs1".into()),
                chr: "chr1".into(),
                start: 90,
                end: 110,
                strand: None,
                alt_info: None,
            },
            VcfRecord {
                chrom: "chr1".into(),
                pos: 100,
                id: Some("rs1".into()),
                ref_allele: "A".into(),
                alt_allele: "G".into(),
            },
        )];
        let expanded = expand_vcf_alt_seq(results, false);
        assert_eq!(expanded.len(), 1);
        assert!(expanded[0].0.alt_info.is_some());
        let info = expanded[0].0.alt_info.as_ref().unwrap();
        assert_eq!(info.ref_allele, "A");
        assert_eq!(info.alt_allele, "G");
        assert_eq!(info.variant_pos, 100);
    }

    #[test]
    fn expand_alt_seq_multi_allelic() {
        let results = vec![(
            Region {
                name: Some("rs1".into()),
                chr: "chr1".into(),
                start: 90,
                end: 110,
                strand: None,
                alt_info: None,
            },
            VcfRecord {
                chrom: "chr1".into(),
                pos: 100,
                id: Some("rs1".into()),
                ref_allele: "A".into(),
                alt_allele: "G,T".into(),
            },
        )];
        let expanded = expand_vcf_alt_seq(results, false);
        assert_eq!(
            expanded.len(),
            2,
            "multi-allelic should produce two entries"
        );
        assert_eq!(expanded[0].0.alt_info.as_ref().unwrap().alt_allele, "G");
        assert_eq!(expanded[1].0.alt_info.as_ref().unwrap().alt_allele, "T");
    }

    #[test]
    fn expand_alt_seq_both_mode() {
        let results = vec![(
            Region {
                name: Some("rs1".into()),
                chr: "chr1".into(),
                start: 90,
                end: 110,
                strand: None,
                alt_info: None,
            },
            VcfRecord {
                chrom: "chr1".into(),
                pos: 100,
                id: Some("rs1".into()),
                ref_allele: "A".into(),
                alt_allele: "G".into(),
            },
        )];
        let expanded = expand_vcf_alt_seq(results, true);
        assert_eq!(expanded.len(), 2, "both mode should produce ref + alt");
        assert!(expanded[0].0.alt_info.is_none(), "first entry is ref");
        assert!(expanded[1].0.alt_info.is_some(), "second entry is alt");
    }
}
