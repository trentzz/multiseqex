use anyhow::{Result, anyhow};
use std::collections::HashMap;
use std::fs::File;
use std::io::{Read, Seek, SeekFrom};

use crate::fai::FaiRecord;
use crate::region::Region;

const BULK_READ_GAP_THRESHOLD: u64 = 8192;

fn region_byte_range(rec: &FaiRecord, r: &Region) -> (u64, u64) {
    let lb = rec.line_bases;
    let lby = rec.line_bytes;
    let start_line = (r.start - 1) / lb;
    let start_col = (r.start - 1) % lb;
    let byte_start = rec.offset + start_line * lby + start_col;
    let end_line = (r.end - 1) / lb;
    let end_col = (r.end - 1) % lb;
    let byte_end = rec.offset + end_line * lby + end_col;
    (byte_start, byte_end)
}

fn strip_and_uppercase(buf: &[u8], expected_bases: usize) -> Vec<u8> {
    let mut seq = Vec::<u8>::with_capacity(expected_bases);
    for &b in buf {
        if b != b'\n' && b != b'\r' {
            seq.push(b.to_ascii_uppercase());
        }
    }
    seq
}

pub fn extract_region(
    f: &mut File,
    fai: &HashMap<String, FaiRecord>,
    r: &Region,
) -> Result<String> {
    if r.start == 0 {
        return Err(anyhow!(
            "region start must be >= 1 (1-based coordinates), got 0 for '{}'",
            r.chr
        ));
    }
    let rec = fai
        .get(&r.chr)
        .ok_or_else(|| anyhow!("Contig '{}' not in index", r.chr))?;
    if rec.line_bases == 0 {
        return Err(anyhow!(
            "Malformed FAI: line_bases is 0 for contig '{}' (would cause division by zero)",
            r.chr
        ));
    }
    let (byte_start, byte_end) = region_byte_range(rec, r);
    let read_len = (byte_end - byte_start + 1) as usize;
    f.seek(SeekFrom::Start(byte_start))?;
    let mut buf = vec![0u8; read_len];
    f.read_exact(&mut buf)?;
    let expected_bases = (r.end - r.start + 1) as usize;
    let seq = strip_and_uppercase(&buf, expected_bases);
    Ok(String::from_utf8(seq)?)
}

pub struct BulkGroup {
    pub indices: Vec<usize>,
    byte_ranges: Vec<(u64, u64)>,
    group_start: u64,
    group_end: u64,
}

pub fn build_bulk_groups(
    regions: &[Region],
    fai: &HashMap<String, FaiRecord>,
) -> Result<Vec<BulkGroup>> {
    let mut indexed: Vec<(usize, u64, u64)> = Vec::with_capacity(regions.len());
    for (i, r) in regions.iter().enumerate() {
        if r.start == 0 {
            return Err(anyhow!(
                "region start must be >= 1 (1-based coordinates), got 0 for '{}'",
                r.chr
            ));
        }
        let rec = fai
            .get(&r.chr)
            .ok_or_else(|| anyhow!("Contig '{}' not in index", r.chr))?;
        if rec.line_bases == 0 {
            return Err(anyhow!(
                "Malformed FAI: line_bases is 0 for contig '{}' (would cause division by zero)",
                r.chr
            ));
        }
        let (bs, be) = region_byte_range(rec, r);
        indexed.push((i, bs, be));
    }
    indexed.sort_by(|a, b| regions[a.0].chr.cmp(&regions[b.0].chr).then(a.1.cmp(&b.1)));
    let mut groups: Vec<BulkGroup> = Vec::new();
    for &(orig_idx, bs, be) in &indexed {
        let contig = regions[orig_idx].chr.as_str();
        let can_merge = groups.last().is_some_and(|g| {
            let last_contig = regions[*g.indices.last().unwrap()].chr.as_str();
            last_contig == contig && bs.saturating_sub(g.group_end) <= BULK_READ_GAP_THRESHOLD
        });
        if can_merge {
            let g = groups.last_mut().unwrap();
            g.indices.push(orig_idx);
            g.byte_ranges.push((bs, be));
            if be > g.group_end {
                g.group_end = be;
            }
        } else {
            groups.push(BulkGroup {
                indices: vec![orig_idx],
                byte_ranges: vec![(bs, be)],
                group_start: bs,
                group_end: be,
            });
        }
    }
    Ok(groups)
}

pub fn extract_bulk_group(
    f: &mut File,
    group: &BulkGroup,
    regions: &[Region],
) -> Result<Vec<(usize, String)>> {
    let total_len = (group.group_end - group.group_start + 1) as usize;
    f.seek(SeekFrom::Start(group.group_start))?;
    let mut buf = vec![0u8; total_len];
    f.read_exact(&mut buf)?;
    let mut results = Vec::with_capacity(group.indices.len());
    for (i, &orig_idx) in group.indices.iter().enumerate() {
        let (bs, be) = group.byte_ranges[i];
        let local_start = (bs - group.group_start) as usize;
        let local_end = (be - group.group_start + 1) as usize;
        let slice = &buf[local_start..local_end];
        let r = &regions[orig_idx];
        let expected_bases = (r.end - r.start + 1) as usize;
        let seq = strip_and_uppercase(slice, expected_bases);
        results.push((orig_idx, String::from_utf8(seq)?));
    }
    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn extract_region_zero_start_errors() {
        let mut tmp = tempfile::NamedTempFile::new().unwrap();
        write!(tmp, ">chr1\nACGTACGT\n").unwrap();
        tmp.flush().unwrap();
        let mut fai = HashMap::new();
        fai.insert(
            "chr1".to_string(),
            FaiRecord {
                length: 8,
                offset: 6,
                line_bases: 8,
                line_bytes: 9,
            },
        );
        let r = Region {
            name: None,
            chr: "chr1".to_string(),
            start: 0,
            end: 5,
            strand: None,
            alt_info: None,
        };
        let mut f = File::open(tmp.path()).unwrap();
        let err = extract_region(&mut f, &fai, &r).unwrap_err();
        assert!(err.to_string().contains("1-based coordinates"));
    }

    #[test]
    fn bulk_read_matches_individual_extraction() {
        let mut tmp = tempfile::NamedTempFile::new().unwrap();
        write!(tmp, ">chr1\nAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCC\nGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTT\n").unwrap();
        tmp.flush().unwrap();
        let mut fai = HashMap::new();
        fai.insert(
            "chr1".to_string(),
            FaiRecord {
                length: 80,
                offset: 6,
                line_bases: 40,
                line_bytes: 41,
            },
        );
        let regions = vec![
            Region {
                name: None,
                chr: "chr1".to_string(),
                start: 1,
                end: 20,
                strand: None,
                alt_info: None,
            },
            Region {
                name: None,
                chr: "chr1".to_string(),
                start: 15,
                end: 40,
                strand: None,
                alt_info: None,
            },
            Region {
                name: None,
                chr: "chr1".to_string(),
                start: 50,
                end: 80,
                strand: None,
                alt_info: None,
            },
        ];
        let mut individual = Vec::new();
        for r in &regions {
            let mut f = File::open(tmp.path()).unwrap();
            individual.push(extract_region(&mut f, &fai, r).unwrap());
        }
        let groups = build_bulk_groups(&regions, &fai).unwrap();
        let mut bulk_results: Vec<(usize, String)> = Vec::new();
        for group in &groups {
            let mut f = File::open(tmp.path()).unwrap();
            bulk_results.extend(extract_bulk_group(&mut f, group, &regions).unwrap());
        }
        bulk_results.sort_by_key(|(idx, _)| *idx);
        assert_eq!(bulk_results.len(), regions.len());
        for (i, (idx, seq)) in bulk_results.iter().enumerate() {
            assert_eq!(*idx, i);
            assert_eq!(
                seq, &individual[i],
                "bulk differs from individual for region {i}"
            );
        }
    }
}

/// Apply an alternate-allele substitution to an extracted reference sequence.
///
/// Given the extracted sequence for a region `[region_start..region_end]`,
/// replaces the REF allele at `variant_pos` with the ALT allele. The
/// variant position is 1-based genomic coordinates. The replacement is
/// case-insensitive on the REF match (the extracted sequence is typically
/// uppercased).
///
/// Returns `Err` if the REF allele at the expected offset does not match
/// the extracted sequence.
pub fn apply_alt_substitution(
    seq: &str,
    region_start: u64,
    info: &crate::region::AltInfo,
) -> Result<String> {
    let offset = (info.variant_pos - region_start) as usize;
    let ref_len = info.ref_allele.len();

    if offset + ref_len > seq.len() {
        return Err(anyhow!(
            "Alt-seq substitution out of bounds: variant at {} with REF len {} \
             exceeds extracted sequence length {} (region starts at {})",
            info.variant_pos,
            ref_len,
            seq.len(),
            region_start
        ));
    }

    let existing = &seq[offset..offset + ref_len];
    if !existing.eq_ignore_ascii_case(&info.ref_allele) {
        return Err(anyhow!(
            "Alt-seq REF mismatch at position {}: expected '{}', found '{}' in extracted sequence",
            info.variant_pos,
            info.ref_allele,
            existing
        ));
    }

    let mut result = String::with_capacity(seq.len() - ref_len + info.alt_allele.len());
    result.push_str(&seq[..offset]);
    result.push_str(&info.alt_allele.to_uppercase());
    result.push_str(&seq[offset + ref_len..]);
    Ok(result)
}

/// Reverse complement a DNA sequence, supporting all IUPAC ambiguity codes.
///
/// Iterates the input as bytes in reverse order, complements each base
/// according to IUPAC rules, and collects to a `Vec<u8>` before converting
/// back to a `String`. This avoids per-character UTF-8 decoding overhead.
/// Unrecognised bytes are left unchanged.
pub fn reverse_complement(seq: &str) -> String {
    let out: Vec<u8> = seq
        .as_bytes()
        .iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            b'R' => b'Y',
            b'Y' => b'R',
            b'S' => b'S',
            b'W' => b'W',
            b'K' => b'M',
            b'M' => b'K',
            b'B' => b'V',
            b'V' => b'B',
            b'D' => b'H',
            b'H' => b'D',
            b'N' => b'N',
            b'a' => b't',
            b't' => b'a',
            b'c' => b'g',
            b'g' => b'c',
            b'r' => b'y',
            b'y' => b'r',
            b's' => b's',
            b'w' => b'w',
            b'k' => b'm',
            b'm' => b'k',
            b'b' => b'v',
            b'v' => b'b',
            b'd' => b'h',
            b'h' => b'd',
            b'n' => b'n',
            other => other,
        })
        .collect();
    // Safety: the input is valid UTF-8 and our mapping preserves ASCII,
    // so the output is always valid UTF-8.
    String::from_utf8(out).expect("reverse_complement produced invalid UTF-8")
}

#[cfg(test)]
mod alt_seq_tests {
    use super::apply_alt_substitution;
    use crate::region::AltInfo;

    fn make_info(ref_allele: &str, alt_allele: &str, variant_pos: u64) -> AltInfo {
        AltInfo {
            ref_allele: ref_allele.to_string(),
            alt_allele: alt_allele.to_string(),
            variant_pos,
        }
    }

    #[test]
    fn snp_substitution() {
        // Region 2-8, variant at pos 5: AACCCGG -> AACGCGG
        let info = make_info("C", "G", 5);
        let result = apply_alt_substitution("AACCCGG", 2, &info).unwrap();
        assert_eq!(result, "AACGCGG");
    }

    #[test]
    fn insertion() {
        // Region 2-8, variant at pos 5: REF=C, ALT=CTG
        // AACCCGG -> AAC[CTG]CGG (offset 3, replace 1 char with 3)
        let info = make_info("C", "CTG", 5);
        let result = apply_alt_substitution("AACCCGG", 2, &info).unwrap();
        assert_eq!(result, "AACCTGCGG");
    }

    #[test]
    fn deletion() {
        // Region 2-10, variant at pos 5: REF=CCG, ALT=C
        // AACCCGGGT -> AAC[C]GGT (offset 3, replace 3 chars with 1)
        let info = make_info("CCG", "C", 5);
        let result = apply_alt_substitution("AACCCGGGT", 2, &info).unwrap();
        assert_eq!(result, "AACCGGT");
    }

    #[test]
    fn ref_mismatch_errors() {
        let info = make_info("T", "G", 5);
        let err = apply_alt_substitution("AACCCGG", 2, &info).unwrap_err();
        assert!(err.to_string().contains("REF mismatch"));
    }

    #[test]
    fn out_of_bounds_errors() {
        let info = make_info("CCCCC", "G", 5);
        let err = apply_alt_substitution("AACCCGG", 2, &info).unwrap_err();
        assert!(err.to_string().contains("out of bounds"));
    }

    #[test]
    fn case_insensitive_ref_match() {
        // Extracted sequence is uppercase, REF in VCF might be lowercase.
        let info = make_info("c", "G", 5);
        let result = apply_alt_substitution("AACCCGG", 2, &info).unwrap();
        assert_eq!(result, "AACGCGG");
    }
}

#[cfg(test)]
mod rc_tests {
    use super::reverse_complement;

    #[test]
    fn rc_standard_bases() {
        assert_eq!(reverse_complement("ACGT"), "ACGT");
        assert_eq!(reverse_complement("AAAA"), "TTTT");
        assert_eq!(reverse_complement("CCCC"), "GGGG");
        assert_eq!(reverse_complement("ATCG"), "CGAT");
    }

    #[test]
    fn rc_iupac_ambiguity_codes() {
        assert_eq!(reverse_complement("R"), "Y");
        assert_eq!(reverse_complement("Y"), "R");
        assert_eq!(reverse_complement("S"), "S");
        assert_eq!(reverse_complement("W"), "W");
        assert_eq!(reverse_complement("K"), "M");
        assert_eq!(reverse_complement("M"), "K");
        assert_eq!(reverse_complement("B"), "V");
        assert_eq!(reverse_complement("V"), "B");
        assert_eq!(reverse_complement("D"), "H");
        assert_eq!(reverse_complement("H"), "D");
        assert_eq!(reverse_complement("N"), "N");
    }

    #[test]
    fn rc_all_iupac_together() {
        assert_eq!(reverse_complement("ACGTRYWSKMBVDHN"), "NDHBVKMSWRYACGT");
    }

    #[test]
    fn rc_empty() {
        assert_eq!(reverse_complement(""), "");
    }

    #[test]
    fn rc_lowercase() {
        assert_eq!(reverse_complement("acgt"), "acgt");
    }
}
