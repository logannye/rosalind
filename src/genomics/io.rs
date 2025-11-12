use anyhow::Result;
use rust_htslib::bam::{self, header::Header, header::HeaderRecord, Writer};
use std::path::Path;

/// Create a BAM writer with a minimal header for a single-reference alignment.
///
/// The caller is responsible for writing alignment records using the returned writer.
pub fn create_bam_writer<P: AsRef<Path>>(
    output_path: P,
    reference_name: &str,
    reference_length: usize,
) -> Result<Writer> {
    let mut header = Header::new();

    let mut hd = HeaderRecord::new(b"HD");
    hd.push_tag(b"VN", &"1.6");
    hd.push_tag(b"SO", &"unknown");
    header.push_record(&hd);

    let mut sq = HeaderRecord::new(b"SQ");
    sq.push_tag(b"SN", reference_name);
    sq.push_tag(b"LN", &(reference_length as i64));
    header.push_record(&sq);

    let writer = bam::Writer::from_path(output_path, &header, bam::Format::Bam)?;
    Ok(writer)
}
