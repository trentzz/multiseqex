//! Multi-sequence extractor for FASTA files using FAI indexing.
pub mod extract;
pub mod fai;
pub mod output;
pub mod region;
pub mod table;
pub mod validate;

// Convenience re-exports for library consumers.
pub use extract::{extract_region, reverse_complement};
pub use fai::{FaiRecord, build_fai, read_fai};
pub use output::wrap_fasta;
pub use region::{Region, deduplicate_regions, parse_region_str, parse_regions_bed, sort_regions};
