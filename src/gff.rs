//! Minimal GFF3/GTF parser for extracting regions from annotation files.
//!
//! GFF3 lines are tab-separated with 9 fields: seqid, source, type, start,
//! end, score, strand, phase, attributes. Coordinates are 1-based inclusive.
//! Lines beginning with `#` are skipped. Features are filtered by type
//! (default: "gene").

use anyhow::{Context, Result, anyhow};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use crate::region::Region;

/// Extract the Name from a GFF3 attributes string (key=value pairs separated by `;`).
/// Falls back to parsing GTF-style attributes (`gene_name "value"`).
fn parse_name_from_attributes(attrs: &str) -> Option<String> {
    // GFF3 style: Name=value
    for part in attrs.split(';') {
        let part = part.trim();
        if let Some(val) = part.strip_prefix("Name=") {
            let val = val.trim();
            if !val.is_empty() {
                return Some(val.to_string());
            }
        }
    }
    // GTF style: gene_name "value"
    for part in attrs.split(';') {
        let part = part.trim();
        if let Some(rest) = part.strip_prefix("gene_name") {
            let rest = rest.trim();
            if let Some(val) = rest.strip_prefix('"') {
                if let Some(val) = val.strip_suffix('"') {
                    if !val.is_empty() {
                        return Some(val.to_string());
                    }
                }
            }
        }
    }
    None
}

/// Parse a strand character from a GFF strand field.
fn parse_strand(s: &str) -> Option<char> {
    match s.trim() {
        "+" => Some('+'),
        "-" => Some('-'),
        "." => Some('.'),
        _ => None,
    }
}

/// Parse a GFF3/GTF file and return regions filtered by feature type.
///
/// GFF3 coordinates are 1-based inclusive. No coordinate conversion is needed.
/// The Name attribute is extracted from the attributes column (column 9).
/// Flanking is not applied here; callers apply it separately if needed.
pub fn parse_regions_gff(path: &Path, feature_type: &str) -> Result<Vec<Region>> {
    let f =
        File::open(path).with_context(|| format!("Cannot open GFF file: {}", path.display()))?;
    let reader = BufReader::new(f);
    let mut regions = Vec::new();

    for (line_num, line_result) in reader.lines().enumerate() {
        let line = line_result
            .with_context(|| format!("I/O error reading GFF file: {}", path.display()))?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = trimmed.split('\t').collect();
        if fields.len() < 9 {
            return Err(anyhow!(
                "GFF line {} has fewer than 9 fields: {}",
                line_num + 1,
                trimmed
            ));
        }

        let ftype = fields[2];
        if ftype != feature_type {
            continue;
        }

        let seqid = fields[0].to_string();
        let start: u64 = fields[3]
            .parse()
            .with_context(|| format!("Bad start at GFF line {}: {}", line_num + 1, fields[3]))?;
        let end: u64 = fields[4]
            .parse()
            .with_context(|| format!("Bad end at GFF line {}: {}", line_num + 1, fields[4]))?;

        if start == 0 {
            return Err(anyhow!(
                "GFF line {} has start=0; GFF uses 1-based coordinates",
                line_num + 1
            ));
        }

        let strand = parse_strand(fields[6]);
        let name = parse_name_from_attributes(fields[8]);

        regions.push(Region {
            name,
            chr: seqid,
            start,
            end,
            strand,
        });
    }
    Ok(regions)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_gff3_gene() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        std::io::Write::write_all(
            &mut tmp.as_file(),
            b"##gff-version 3\nchr1\t.\tgene\t100\t200\t.\t+\t.\tID=gene1;Name=TP53\n\
              chr1\t.\texon\t100\t150\t.\t+\t.\tParent=gene1\n\
              chr2\t.\tgene\t300\t400\t.\t-\t.\tID=gene2;Name=BRCA1\n",
        )
        .unwrap();
        let regions = parse_regions_gff(tmp.path(), "gene").unwrap();
        assert_eq!(regions.len(), 2);
        assert_eq!(regions[0].chr, "chr1");
        assert_eq!(regions[0].start, 100);
        assert_eq!(regions[0].end, 200);
        assert_eq!(regions[0].strand, Some('+'));
        assert_eq!(regions[0].name.as_deref(), Some("TP53"));
        assert_eq!(regions[1].chr, "chr2");
        assert_eq!(regions[1].start, 300);
        assert_eq!(regions[1].end, 400);
        assert_eq!(regions[1].strand, Some('-'));
        assert_eq!(regions[1].name.as_deref(), Some("BRCA1"));
    }

    #[test]
    fn parse_gff3_exon_filter() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        std::io::Write::write_all(
            &mut tmp.as_file(),
            b"chr1\t.\tgene\t100\t200\t.\t+\t.\tName=TP53\n\
              chr1\t.\texon\t100\t150\t.\t+\t.\tName=exon1\n\
              chr1\t.\texon\t160\t200\t.\t+\t.\tName=exon2\n",
        )
        .unwrap();
        let regions = parse_regions_gff(tmp.path(), "exon").unwrap();
        assert_eq!(regions.len(), 2);
        assert_eq!(regions[0].name.as_deref(), Some("exon1"));
        assert_eq!(regions[1].name.as_deref(), Some("exon2"));
    }

    #[test]
    fn parse_gtf_gene_name() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        std::io::Write::write_all(
            &mut tmp.as_file(),
            b"chr1\tensembl\tgene\t100\t200\t.\t+\t.\tgene_id \"ENSG001\"; gene_name \"TP53\"\n",
        )
        .unwrap();
        let regions = parse_regions_gff(tmp.path(), "gene").unwrap();
        assert_eq!(regions.len(), 1);
        assert_eq!(regions[0].name.as_deref(), Some("TP53"));
    }

    #[test]
    fn skip_comments_and_blanks() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        std::io::Write::write_all(
            &mut tmp.as_file(),
            b"##gff-version 3\n# comment\n\nchr1\t.\tgene\t1\t10\t.\t.\t.\tName=g1\n",
        )
        .unwrap();
        let regions = parse_regions_gff(tmp.path(), "gene").unwrap();
        assert_eq!(regions.len(), 1);
    }

    #[test]
    fn dot_strand() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        std::io::Write::write_all(
            &mut tmp.as_file(),
            b"chr1\t.\tgene\t1\t10\t.\t.\t.\tName=g1\n",
        )
        .unwrap();
        let regions = parse_regions_gff(tmp.path(), "gene").unwrap();
        assert_eq!(regions[0].strand, Some('.'));
    }

    #[test]
    fn no_name_gives_none() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        std::io::Write::write_all(
            &mut tmp.as_file(),
            b"chr1\t.\tgene\t1\t10\t.\t+\t.\tID=gene1\n",
        )
        .unwrap();
        let regions = parse_regions_gff(tmp.path(), "gene").unwrap();
        assert!(regions[0].name.is_none());
    }

    #[test]
    fn start_zero_errors() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        std::io::Write::write_all(
            &mut tmp.as_file(),
            b"chr1\t.\tgene\t0\t10\t.\t+\t.\tName=g1\n",
        )
        .unwrap();
        let err = parse_regions_gff(tmp.path(), "gene").unwrap_err();
        assert!(err.to_string().contains("1-based"));
    }

    #[test]
    fn fewer_than_9_fields_errors() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        std::io::Write::write_all(&mut tmp.as_file(), b"chr1\t.\tgene\t1\t10\t.\t+\n").unwrap();
        let err = parse_regions_gff(tmp.path(), "gene").unwrap_err();
        assert!(err.to_string().contains("fewer than 9 fields"));
    }

    #[test]
    fn parse_name_gff3_style() {
        assert_eq!(
            parse_name_from_attributes("ID=gene1;Name=TP53"),
            Some("TP53".into())
        );
    }

    #[test]
    fn parse_name_gtf_style() {
        assert_eq!(
            parse_name_from_attributes("gene_id \"ENSG001\"; gene_name \"TP53\""),
            Some("TP53".into())
        );
    }

    #[test]
    fn parse_name_none() {
        assert_eq!(
            parse_name_from_attributes("ID=gene1;Dbxref=GeneID:123"),
            None
        );
    }
}
