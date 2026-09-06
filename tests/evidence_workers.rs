use std::path::PathBuf;
use std::process::Command;

use rosalind::evidence::CANONICAL_TILE_BASES;
use rosalind::provenance::RunManifest;
use rust_htslib::bam::{
    self,
    header::HeaderRecord,
    record::{Cigar, CigarString},
    Header, Record,
};

#[test]
fn eight_canonical_partitions_are_identical_across_worker_and_budget_matrix() {
    let root = std::env::temp_dir().join(format!(
        "rosalind-worker-matrix-{}-{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    std::fs::create_dir(&root).unwrap();
    struct Cleanup(PathBuf);
    impl Drop for Cleanup {
        fn drop(&mut self) {
            let _ = std::fs::remove_dir_all(&self.0);
        }
    }
    let _cleanup = Cleanup(root.clone());
    let length = 8 * CANONICAL_TILE_BASES;
    let fasta = root.join("reference.fa");
    let bam = root.join("reads.bam");
    let bed = root.join("sites.bed");
    std::fs::write(&fasta, format!(">chr1\n{}\n", "A".repeat(length as usize))).unwrap();
    std::fs::write(
        root.join("reference.fa.fai"),
        format!("chr1\t{length}\t6\t{length}\t{}\n", length + 1),
    )
    .unwrap();
    let mut header = Header::new();
    header.push_record(
        HeaderRecord::new(b"HD")
            .push_tag(b"VN", "1.6")
            .push_tag(b"SO", "coordinate"),
    );
    header.push_record(
        HeaderRecord::new(b"SQ")
            .push_tag(b"SN", "chr1")
            .push_tag(b"LN", length),
    );
    let mut writer = bam::Writer::from_path(&bam, &header, bam::Format::Bam).unwrap();
    let mut regions = String::new();
    for tile in 0..8 {
        let position = tile * CANONICAL_TILE_BASES + 3;
        regions.push_str(&format!(
            "chr1\t{position}\t{}\tlocus-{tile}\n",
            position + 1
        ));
        for allele in 0..4 {
            let mut record = Record::new();
            record.set(
                format!("r-{tile}-{allele}").as_bytes(),
                Some(&CigarString(vec![Cigar::Match(8)])),
                &[b"ACGT"[allele]; 8],
                &[30; 8],
            );
            record.set_tid(0);
            record.set_pos(i64::from(position - 3));
            record.set_flags(if allele % 2 == 0 { 0 } else { 16 });
            record.set_mapq(60);
            writer.write(&record).unwrap();
        }
    }
    drop(writer);
    bam::index::build(&bam, None::<&PathBuf>, bam::index::Type::Bai, 1).unwrap();
    std::fs::write(&bed, regions).unwrap();
    let mut canonical_bytes = None;
    let mut canonical_identity = None;
    for (workers, budget) in [(1, 256), (2, 384), (8, 1024)] {
        let output = root.join(format!("workers-{workers}.arrow"));
        let result = Command::new(env!("CARGO_BIN_EXE_rosalind"))
            .args(["analyze", "evidence", "--alignments"])
            .arg(&bam)
            .arg("--reference")
            .arg(&fasta)
            .arg("--regions")
            .arg(&bed)
            .args([
                "--format",
                "arrow-ipc",
                "--workers",
                &workers.to_string(),
                "--memory-budget-mb",
                &budget.to_string(),
            ])
            .arg("--cache-dir")
            .arg(root.join(format!("cache-{workers}")))
            .arg("-o")
            .arg(&output)
            .output()
            .unwrap();
        assert!(
            result.status.success(),
            "workers={workers}, budget={budget}, exit={:?}: {}",
            result.status.code(),
            String::from_utf8_lossy(&result.stderr)
        );
        let receipt = RunManifest::from_canonical_json(
            &std::fs::read_to_string(format!("{}.manifest.json", output.display())).unwrap(),
        )
        .unwrap();
        assert_eq!(receipt.self_hash_ok(), Some(true));
        assert_eq!(receipt.params["run_status"], "completed");
        assert_eq!(
            receipt.measurements["execution.worker_count"],
            workers.to_string()
        );
        assert_eq!(receipt.measurements["execution.computed_partitions"], "8");
        assert_eq!(receipt.measurements["execution.emitted_loci"], "8");
        let identity = receipt.params["evidence.science_blake3"].clone();
        let bytes = std::fs::read(&output).unwrap();
        if let Some(canonical) = &canonical_bytes {
            assert_eq!(&bytes, canonical);
        } else {
            canonical_bytes = Some(bytes);
        }
        if let Some(canonical) = &canonical_identity {
            assert_eq!(&identity, canonical);
        } else {
            canonical_identity = Some(identity);
        }
    }
}
