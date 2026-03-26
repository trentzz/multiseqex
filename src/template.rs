//! Name template expansion for FASTA/FASTQ headers.
//!
//! Supported placeholders: `{chr}`, `{start}`, `{end}`, `{name}`, `{length}`,
//! `{index}`, `{strand}`.

use crate::region::Region;

/// Expand a name template string, replacing placeholders with region values.
///
/// `index` is the 1-based region index in the output.
pub fn expand_template(template: &str, region: &Region, index: usize) -> String {
    let length = region.end.saturating_sub(region.start) + 1;
    let strand = match region.strand {
        Some(c) => String::from(c),
        None => ".".to_string(),
    };
    let name = region.name.as_deref().unwrap_or("");

    template
        .replace("{chr}", &region.chr)
        .replace("{start}", &region.start.to_string())
        .replace("{end}", &region.end.to_string())
        .replace("{name}", name)
        .replace("{length}", &length.to_string())
        .replace("{index}", &index.to_string())
        .replace("{strand}", &strand)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_region(name: Option<&str>, strand: Option<char>) -> Region {
        Region {
            name: name.map(|s| s.to_string()),
            chr: "chr1".into(),
            start: 100,
            end: 200,
            strand,
        }
    }

    #[test]
    fn all_placeholders() {
        let r = make_region(Some("gene1"), Some('+'));
        let result = expand_template(
            "{chr}:{start}-{end} {name} len={length} idx={index} strand={strand}",
            &r,
            5,
        );
        assert_eq!(result, "chr1:100-200 gene1 len=101 idx=5 strand=+");
    }

    #[test]
    fn no_name_gives_empty() {
        let r = make_region(None, None);
        let result = expand_template("{name}", &r, 1);
        assert_eq!(result, "");
    }

    #[test]
    fn no_strand_gives_dot() {
        let r = make_region(None, None);
        let result = expand_template("{strand}", &r, 1);
        assert_eq!(result, ".");
    }

    #[test]
    fn minus_strand() {
        let r = make_region(None, Some('-'));
        let result = expand_template("{strand}", &r, 1);
        assert_eq!(result, "-");
    }

    #[test]
    fn plain_text_no_placeholders() {
        let r = make_region(None, None);
        let result = expand_template("hello world", &r, 1);
        assert_eq!(result, "hello world");
    }

    #[test]
    fn index_is_one_based() {
        let r = make_region(None, None);
        let result = expand_template("{index}", &r, 1);
        assert_eq!(result, "1");
    }

    #[test]
    fn length_calculation() {
        let r = Region {
            name: None,
            chr: "chr1".into(),
            start: 1,
            end: 10,
            strand: None,
        };
        let result = expand_template("{length}", &r, 1);
        assert_eq!(result, "10");
    }
}
