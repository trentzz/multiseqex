//! Multi-sequence extractor for FASTA files using FAI indexing.
pub mod extract;
pub mod fai;
pub mod gff;
pub mod intervals;
pub mod mask;
pub mod noindex;
pub mod output;
pub mod region;
pub mod stats;
pub mod table;
pub mod template;
pub mod transform;
pub mod validate;
pub mod vcf;

// Convenience re-exports for library consumers.
pub use extract::{extract_region, reverse_complement};
pub use fai::{FaiRecord, build_fai, read_fai};
pub use gff::parse_regions_gff;
pub use intervals::{intersect_regions, subtract_regions};
pub use mask::{MaskIndex, MaskMode};
pub use output::wrap_fasta;
pub use region::{
    Region, deduplicate_regions, merge_regions, parse_region_str, parse_regions_bed,
    resolve_flanks, sort_regions, tile_regions,
};
pub use stats::{RegionStats, compute_stats};
pub use template::expand_template;
pub use transform::TransformConfig;
pub use validate::{is_gzip, resolve_bgzip};
pub use vcf::{VcfRecord, parse_regions_vcf};
