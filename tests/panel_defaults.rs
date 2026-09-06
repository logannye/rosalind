use rosalind::evidence::{
    EvidenceEngine, EvidenceRequest, PanelQcAnalyzer, PanelTarget, DEFAULT_MIN_CALLABLE_DEPTH,
};
use rust_htslib::bam::record::{Cigar, CigarString, Record};
use rust_htslib::bam::{self, header::HeaderRecord, Header};
use std::path::PathBuf;
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

static FIXTURE_ID: AtomicU64 = AtomicU64::new(0);
const DEPTHS: [u64; 4] = [9, 10, 19, 20];

struct Fixture {
    dir: PathBuf,
    reference: PathBuf,
    bam: PathBuf,
    bed: PathBuf,
}

impl Fixture {
    fn new() -> Self {
        let dir = std::env::temp_dir().join(format!(
            "rosalind-panel-defaults-{}-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos(),
            FIXTURE_ID.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir(&dir).unwrap();
        let reference = dir.join("reference.fa");
        std::fs::write(&reference, ">chr1\nAAAA\n").unwrap();
        std::fs::write(dir.join("reference.fa.fai"), "chr1\t4\t6\t4\t5\n").unwrap();
        let bam = dir.join("reads.bam");
        let mut header = Header::new();
        header.push_record(
            HeaderRecord::new(b"HD")
                .push_tag(b"VN", "1.6")
                .push_tag(b"SO", "coordinate"),
        );
        header.push_record(
            HeaderRecord::new(b"SQ")
                .push_tag(b"SN", "chr1")
                .push_tag(b"LN", 4),
        );
        let mut writer = bam::Writer::from_path(&bam, &header, bam::Format::Bam).unwrap();
        for (position, depth) in DEPTHS.into_iter().enumerate() {
            for read in 0..depth {
                let mut record = Record::new();
                record.set(
                    format!("p{position}-r{read}").as_bytes(),
                    Some(&CigarString(vec![Cigar::Match(1)])),
                    b"A",
                    &[30],
                );
                record.set_tid(0);
                record.set_pos(position as i64);
                record.set_mapq(60);
                record.set_flags(0);
                writer.write(&record).unwrap();
            }
        }
        drop(writer);
        bam::index::build(&bam, None::<&PathBuf>, bam::index::Type::Bai, 1).unwrap();
        let bed = dir.join("targets.bed");
        std::fs::write(
            &bed,
            "chr1\t0\t1\tdepth9\nchr1\t1\t2\tdepth10\nchr1\t2\t3\tdepth19\nchr1\t3\t4\tdepth20\nchr1\t0\t4\tall\n",
        )
        .unwrap();
        Self {
            dir,
            reference,
            bam,
            bed,
        }
    }

    fn rust_panel(&self, threshold: Option<u64>) -> (PanelQcAnalyzer, Vec<u8>) {
        let mut engine =
            EvidenceEngine::open(EvidenceRequest::new(&self.bam, &self.reference)).unwrap();
        let targets = PanelTarget::from_bed(&self.bed, engine.contigs()).unwrap();
        let mut panel = PanelQcAnalyzer::new(targets).unwrap();
        if let Some(threshold) = threshold {
            panel = panel.with_min_callable_depth(threshold);
        }
        engine.set_selection(panel.selection()).unwrap();
        engine.run(&mut panel).unwrap();
        let mut output = Vec::new();
        panel.write_tsv(&mut output, engine.contigs()).unwrap();
        (panel, output)
    }

    fn cli_panel(&self, threshold: Option<u64>) -> Vec<u8> {
        let mut command = Command::new(env!("CARGO_BIN_EXE_rosalind"));
        command
            .args(["analyze", "panel-qc", "--alignments"])
            .arg(&self.bam)
            .arg("--reference")
            .arg(&self.reference)
            .arg("--regions")
            .arg(&self.bed);
        if let Some(threshold) = threshold {
            command
                .arg("--min-callable-depth")
                .arg(threshold.to_string());
        }
        let output = command.output().unwrap();
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        output.stdout
    }
}

impl Drop for Fixture {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.dir);
    }
}

#[test]
fn rust_and_cli_share_default_and_explicit_panel_callability_thresholds() {
    let fixture = Fixture::new();
    assert_eq!(DEFAULT_MIN_CALLABLE_DEPTH, 10);
    for (threshold, expected_threshold, expected_callable) in
        [(None, 10, [0, 1, 1, 1, 3]), (Some(20), 20, [0, 0, 0, 1, 1])]
    {
        let (panel, rust_output) = fixture.rust_panel(threshold);
        assert_eq!(panel.min_callable_depth, expected_threshold);
        assert_eq!(
            panel
                .summaries()
                .iter()
                .map(|summary| summary.callable_positions)
                .collect::<Vec<_>>(),
            expected_callable
        );
        for (summary, depth) in panel.summaries().iter().zip(DEPTHS) {
            assert_eq!(summary.callable_depth_sum, depth);
            assert_eq!(summary.base_quality_sum, depth * 30);
            assert_eq!(summary.mapping_quality_sum, depth * 60);
        }
        let whole_panel = &panel.summaries()[4];
        assert_eq!(whole_panel.target_length(), 4);
        assert_eq!(whole_panel.callable_depth_sum, 58);
        assert_eq!(whole_panel.breadth_positions, [4, 3, 1, 0]);
        assert_eq!(fixture.cli_panel(threshold), rust_output);
    }
}
