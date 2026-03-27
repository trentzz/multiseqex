//! Sequence transforms applied after extraction and reverse complement,
//! before output.
//!
//! Transforms: `--to-rna`, `--uppercase`, `--lowercase`, `--translate`.

/// Convert DNA to RNA by replacing T with U (and t with u).
pub fn to_rna(seq: &str) -> String {
    seq.as_bytes()
        .iter()
        .map(|&b| match b {
            b'T' => b'U',
            b't' => b'u',
            other => other,
        })
        .map(|b| b as char)
        .collect()
}

/// Force all bases to uppercase.
pub fn to_uppercase(seq: &str) -> String {
    seq.to_ascii_uppercase()
}

/// Force all bases to lowercase.
pub fn to_lowercase(seq: &str) -> String {
    seq.to_ascii_lowercase()
}

/// Standard genetic code translation table.
/// Translates a three-letter codon (uppercase) to a single amino acid character.
/// Stop codons produce '*'. Unknown codons produce 'X'.
fn translate_codon(codon: &[u8]) -> u8 {
    match codon {
        b"TTT" | b"TTC" => b'F',
        b"TTA" | b"TTG" => b'L',
        b"CTT" | b"CTC" | b"CTA" | b"CTG" => b'L',
        b"ATT" | b"ATC" | b"ATA" => b'I',
        b"ATG" => b'M',
        b"GTT" | b"GTC" | b"GTA" | b"GTG" => b'V',
        b"TCT" | b"TCC" | b"TCA" | b"TCG" => b'S',
        b"CCT" | b"CCC" | b"CCA" | b"CCG" => b'P',
        b"ACT" | b"ACC" | b"ACA" | b"ACG" => b'T',
        b"GCT" | b"GCC" | b"GCA" | b"GCG" => b'A',
        b"TAT" | b"TAC" => b'Y',
        b"TAA" | b"TAG" | b"TGA" => b'*',
        b"CAT" | b"CAC" => b'H',
        b"CAA" | b"CAG" => b'Q',
        b"AAT" | b"AAC" => b'N',
        b"AAA" | b"AAG" => b'K',
        b"GAT" | b"GAC" => b'D',
        b"GAA" | b"GAG" => b'E',
        b"TGT" | b"TGC" => b'C',
        b"TGG" => b'W',
        b"CGT" | b"CGC" | b"CGA" | b"CGG" => b'R',
        b"AGT" | b"AGC" => b'S',
        b"AGA" | b"AGG" => b'R',
        b"GGT" | b"GGC" | b"GGA" | b"GGG" => b'G',
        _ => b'X',
    }
}

/// Translate a DNA sequence to amino acids using the standard genetic code.
///
/// Reading frame starts from position 1 (the first base). Codons are read
/// in non-overlapping triplets. Any trailing bases that do not form a
/// complete codon are ignored. The input is uppercased before translation.
pub fn translate(seq: &str) -> String {
    let upper = seq.to_ascii_uppercase();
    let bytes = upper.as_bytes();
    let mut result = Vec::with_capacity(bytes.len() / 3);
    for chunk in bytes.chunks_exact(3) {
        result.push(translate_codon(chunk));
    }
    String::from_utf8(result).expect("translate produced invalid UTF-8")
}

/// Configuration for sequence transforms.
#[derive(Debug, Clone, Default)]
pub struct TransformConfig {
    pub to_rna: bool,
    pub uppercase: bool,
    pub lowercase: bool,
    pub translate: bool,
}

impl TransformConfig {
    /// Returns true if any transform is enabled.
    pub fn any_active(&self) -> bool {
        self.to_rna || self.uppercase || self.lowercase || self.translate
    }

    /// Apply all configured transforms to a sequence.
    ///
    /// Order: translate or to_rna first, then case transforms.
    pub fn apply(&self, seq: &str) -> String {
        let mut s = seq.to_string();
        if self.translate {
            s = translate(&s);
        } else if self.to_rna {
            s = to_rna(&s);
        }
        if self.uppercase {
            s = to_uppercase(&s);
        } else if self.lowercase {
            s = to_lowercase(&s);
        }
        s
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // ── to_rna ──────────────────────────────────────────────────────────

    #[test]
    fn rna_converts_t_to_u() {
        assert_eq!(to_rna("ACGT"), "ACGU");
    }

    #[test]
    fn rna_preserves_case() {
        assert_eq!(to_rna("acgt"), "acgu");
    }

    #[test]
    fn rna_no_t() {
        assert_eq!(to_rna("ACGN"), "ACGN");
    }

    #[test]
    fn rna_empty() {
        assert_eq!(to_rna(""), "");
    }

    // ── to_uppercase / to_lowercase ─────────────────────────────────────

    #[test]
    fn uppercase_works() {
        assert_eq!(to_uppercase("acgt"), "ACGT");
    }

    #[test]
    fn lowercase_works() {
        assert_eq!(to_lowercase("ACGT"), "acgt");
    }

    // ── translate ───────────────────────────────────────────────────────

    #[test]
    fn translate_atg_start() {
        assert_eq!(translate("ATG"), "M");
    }

    #[test]
    fn translate_stop_codons() {
        assert_eq!(translate("TAA"), "*");
        assert_eq!(translate("TAG"), "*");
        assert_eq!(translate("TGA"), "*");
    }

    #[test]
    fn translate_multiple_codons() {
        // ATG GCT TAA = M A *
        assert_eq!(translate("ATGGCTTAA"), "MA*");
    }

    #[test]
    fn translate_ignores_incomplete_codon() {
        assert_eq!(translate("ATGC"), "M");
    }

    #[test]
    fn translate_lowercase_input() {
        assert_eq!(translate("atg"), "M");
    }

    #[test]
    fn translate_empty() {
        assert_eq!(translate(""), "");
    }

    #[test]
    fn translate_known_protein() {
        // M=ATG, F=TTT, L=CTG, *=TAA
        assert_eq!(translate("ATGTTTCTGTAA"), "MFL*");
    }

    #[test]
    fn translate_unknown_codon() {
        // NNN is not a standard codon, should produce X.
        assert_eq!(translate("NNN"), "X");
    }

    // ── TransformConfig ─────────────────────────────────────────────────

    #[test]
    fn config_no_transforms() {
        let cfg = TransformConfig::default();
        assert!(!cfg.any_active());
        assert_eq!(cfg.apply("ACGT"), "ACGT");
    }

    #[test]
    fn config_rna_and_lowercase() {
        let cfg = TransformConfig {
            to_rna: true,
            lowercase: true,
            ..Default::default()
        };
        assert!(cfg.any_active());
        assert_eq!(cfg.apply("ACGT"), "acgu");
    }

    #[test]
    fn config_translate_and_uppercase() {
        let cfg = TransformConfig {
            translate: true,
            uppercase: true,
            ..Default::default()
        };
        assert_eq!(cfg.apply("atg"), "M");
    }
}
