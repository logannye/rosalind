//! External factory/lifecycle integration tests use real indexed BAM and datasets.
use rosalind::contract::{
    AnalyzerIdentity, EnforcementMode, OutputPolicy, ProducerIdentity, ReplayInvocation,
};
use rosalind::core::cancellation::CancellationToken;
use rosalind::dataset::{
    publish_evidence_dataset, run_dataset_with_snapshot, DatasetOptions, DatasetQuery,
    DatasetReadLimits, DescriptorLimits, VerifiedInputSession,
};
use rosalind::evidence::*;
use rosalind::provenance::{blake3_file, RunManifest};
use rust_htslib::bam::{
    self,
    record::{Cigar, CigarString},
};
use std::{
    collections::BTreeMap,
    fs,
    io::{self, Write},
    path::{Path, PathBuf},
    sync::{
        atomic::{AtomicU64, Ordering},
        Mutex,
    },
};
static SERIAL: Mutex<()> = Mutex::new(());
static NEXT: AtomicU64 = AtomicU64::new(0);
struct Fixture {
    root: PathBuf,
    bam: PathBuf,
    reference: PathBuf,
    bed: PathBuf,
}
impl Fixture {
    fn new() -> Self {
        let root = std::env::temp_dir().join(format!(
            "rosalind-artifact-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir_all(&root).unwrap();
        let reference = root.join("ref.fa");
        fs::write(&reference, format!(">chr1\n{}\n", "A".repeat(40000))).unwrap();
        fs::write(root.join("ref.fa.fai"), "chr1\t40000\t6\t40000\t40001\n").unwrap();
        let bam = root.join("reads.bam");
        let mut header = bam::Header::new();
        header.push_record(
            bam::header::HeaderRecord::new(b"SQ")
                .push_tag(b"SN", "chr1")
                .push_tag(b"LN", 40000),
        );
        let mut writer = bam::Writer::from_path(&bam, &header, bam::Format::Bam).unwrap();
        for (index, pos) in [0, 1020, 16380, 32760].into_iter().enumerate() {
            let mut record = bam::Record::new();
            record.set(
                format!("read{index}").as_bytes(),
                Some(&CigarString(vec![Cigar::Match(20)])),
                b"AACAAAAAACAAAAAAAAAA",
                &[30; 20],
            );
            record.set_tid(0);
            record.set_pos(pos);
            record.set_mapq(60);
            record.set_flags(0);
            writer.write(&record).unwrap();
        }
        drop(writer);
        bam::index::build(&bam, None::<&PathBuf>, bam::index::Type::Bai, 1).unwrap();
        let bed = root.join("selected.bed");
        fs::write(
            &bed,
            "chr1\t0\t2051\nchr1\t16380\t16400\nchr1\t32760\t32780\n",
        )
        .unwrap();
        Self {
            root,
            bam,
            reference,
            bed,
        }
    }
    fn request(&self) -> EvidenceRequest {
        let mut request = EvidenceRequest::new(&self.bam, &self.reference);
        request.fields = fields();
        request.alignment_index = Some(self.root.join("reads.bam.bai"));
        request.reference_fai = Some(self.root.join("ref.fa.fai"));
        request
    }
    fn spec(&self, name: &str) -> EvidenceArtifactSpec {
        EvidenceArtifactSpec::new(
            EvidenceArtifactSource::Native {
                request: self.request(),
                selection: ArtifactSelection::Bed(self.bed.clone()),
            },
            self.root.join(name),
            ProducerIdentity {
                name: "independent-analyzer".into(),
                version: "1.2.3".into(),
                repository: None,
                binary: "independent-analyzer".into(),
            },
            AnalyzerIdentity::new("batch-summary", "2"),
            ReplayInvocation::new(["run"]).option("--scale", 1),
        )
    }
    fn dataset(&self) -> PathBuf {
        let session = VerifiedInputSession::open(vec![
            ("alignments".into(), self.bam.clone()),
            ("alignment-index".into(), self.root.join("reads.bam.bai")),
            ("reference".into(), self.reference.clone()),
            ("reference-fai".into(), self.root.join("ref.fa.fai")),
            ("regions".into(), self.bed.clone()),
        ])
        .unwrap();
        let mut request = self.request();
        request.fields = EvidenceFields::ALL_SUPPORTED;
        let mut engine = EvidenceEngine::open(request).unwrap();
        engine
            .set_selection(EvidenceSelection::from_bed(&self.bed, engine.contigs()).unwrap())
            .unwrap();
        let namespace = session.dataset_namespace(&engine).unwrap();
        let outcome = run_dataset_with_snapshot(
            &mut engine,
            &namespace,
            &DatasetOptions {
                cache_dir: self.root.join("cache"),
                resume: false,
                workers: 1,
            },
            &mut EvidenceCallback::with_fields(
                |_: &EvidenceBatch| Ok(()),
                0,
                EvidenceFields::ALL_SUPPORTED,
            ),
            session.snapshot(),
        )
        .unwrap();
        publish_evidence_dataset(&engine, &outcome, &session, DescriptorLimits::default()).unwrap()
    }
}
impl Drop for Fixture {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.root);
    }
}
fn fields() -> EvidenceFields {
    EvidenceFields::DEPTHS.union(EvidenceFields::ALLELES)
}
#[derive(Clone, Default)]
enum Failure {
    #[default]
    None,
    Create,
    Mismatch,
    Batch,
    Finish,
    Flush,
    Resource,
    Cancel(CancellationToken),
    Mutate(PathBuf),
}
struct Factory {
    failure: Failure,
    unknown: bool,
}
impl Default for Factory {
    fn default() -> Self {
        Self {
            failure: Failure::None,
            unknown: false,
        }
    }
}
impl EvidenceArtifactFactory for Factory {
    fn requirements(&self) -> EvidenceRequirements {
        EvidenceRequirements {
            fields: fields(),
            requires_reference: true,
            context_bases: 0,
            retained_bytes: (!self.unknown).then_some(4096),
        }
    }
    fn params(&self) -> BTreeMap<String, String> {
        BTreeMap::from([
            ("scale".into(), "1".into()),
            ("metric".into(), "canonical-batch-depth".into()),
        ])
    }
    fn create<'a>(
        &'a mut self,
        out: &'a mut dyn Write,
    ) -> Result<Box<dyn EvidenceAnalyzer + 'a>, EvidenceError> {
        writeln!(out, "#contig\tstart\trows\tdepth")?;
        if matches!(self.failure, Failure::Create) {
            return Err(EvidenceError::Analyzer("factory failed".into()));
        }
        let mut requirements = self.requirements();
        if matches!(self.failure, Failure::Mismatch) {
            requirements.retained_bytes = Some(0);
        }
        Ok(Box::new(Consumer {
            out,
            requirements,
            failure: self.failure.clone(),
        }))
    }
}
struct Consumer<'a> {
    out: &'a mut dyn Write,
    requirements: EvidenceRequirements,
    failure: Failure,
}
impl EvidenceAnalyzer for Consumer<'_> {
    fn requirements(&self) -> EvidenceRequirements {
        self.requirements.clone()
    }
    fn on_batch(&mut self, batch: &EvidenceBatch) -> Result<(), EvidenceError> {
        let depth: u64 = batch
            .rows()
            .map(|row| row.depths.unwrap().callable_depth)
            .sum();
        writeln!(
            self.out,
            "{}\t{}\t{}\t{depth}",
            batch.contig,
            batch.row(0).unwrap().position,
            batch.len()
        )?;
        match &self.failure {
            Failure::Batch => return Err(EvidenceError::Analyzer("batch failed".into())),
            Failure::Resource => {
                return Err(EvidenceError::RecordLimit(
                    "injected declared envelope breach".into(),
                ))
            }
            Failure::Cancel(token) => token.cancel(),
            Failure::Mutate(path) => {
                fs::write(path, b"mutated during consumer")?;
            }
            _ => {}
        }
        Ok(())
    }
    fn finish(&mut self) -> Result<(), EvidenceError> {
        match self.failure {
            Failure::Finish => Err(EvidenceError::Analyzer("finish failed".into())),
            Failure::Flush => Err(EvidenceError::Io(io::Error::other("encoder flush failed"))),
            _ => {
                self.out.flush()?;
                Ok(())
            }
        }
    }
}
fn manifest(path: &Path) -> RunManifest {
    RunManifest::from_canonical_json(&fs::read_to_string(path).unwrap()).unwrap()
}

#[test]
fn canonical_batches_and_science_are_identical_across_three_budgets_and_tiles() {
    let _guard = SERIAL.lock().unwrap();
    let f = Fixture::new();
    let mut previous = None;
    let mut science = None;
    for (index, (budget, width)) in [(96, 1), (128, 256), (256, 16384)].into_iter().enumerate() {
        let mut spec = f.spec(&format!("out{index}.tsv"));
        spec.enforcement = EnforcementMode::Cooperative;
        if let EvidenceArtifactSource::Native { request, .. } = &mut spec.source {
            request.execution.memory_budget_bytes = Some(budget << 20);
            request.execution.max_microtile_bases = width;
        }
        let outcome = run_evidence_artifact(&mut Factory::default(), spec).unwrap();
        assert!(outcome.completed);
        assert_eq!(outcome.stats.emitted_loci, 2091);
        let bytes = fs::read(&outcome.output).unwrap();
        if let Some(previous) = &previous {
            assert_eq!(previous, &bytes);
        }
        previous = Some(bytes);
        let receipt = manifest(&outcome.manifest);
        assert_eq!(receipt.self_hash_ok(), Some(true));
        assert_eq!(receipt.measurement_hash_ok(), Some(true));
        assert_eq!(receipt.tool_version, "1.2.3");
        assert_eq!(receipt.params["replay.kind"], "external-analyzer");
        assert_eq!(receipt.params["replay_schema"], "3");
        assert_eq!(receipt.params["contract.canonical_batch_rows"], "1024");
        assert_eq!(
            receipt
                .get_recorded("analyzer.max_additional_bytes")
                .map(String::as_str),
            Some("4096")
        );
        assert_eq!(receipt.params["analyzer.required_fields"], "3");
        assert_eq!(receipt.params["analyzer.requires_reference"], "true");
        assert_eq!(receipt.params["analyzer.context_bases"], "0");
        assert_eq!(
            receipt
                .get_recorded("peak_rss_bytes")
                .unwrap()
                .parse::<u64>()
                .unwrap(),
            outcome.peak_rss_bytes
        );
        assert_eq!(
            receipt
                .get_recorded("memory_budget_bytes")
                .unwrap()
                .parse::<u64>()
                .unwrap(),
            budget << 20
        );
        assert_eq!(
            receipt.outputs[0].blake3,
            blake3_file(&outcome.output).unwrap()
        );
        if let Some(science) = &science {
            assert_eq!(science, &receipt.params["science.blake3"]);
        }
        science = Some(receipt.params["science.blake3"].clone());
    }
    let text = String::from_utf8(previous.unwrap()).unwrap();
    assert!(text.contains("chr1\t0\t1024\t"));
    assert!(text.contains("chr1\t1024\t1024\t"));
    assert!(text.contains("chr1\t2048\t7\t"));
}

#[test]
fn native_and_projected_persisted_source_outputs_and_science_match_after_relocation() {
    let _guard = SERIAL.lock().unwrap();
    let f = Fixture::new();
    let native = run_evidence_artifact(&mut Factory::default(), f.spec("native.tsv")).unwrap();
    let parent = f.dataset();
    let old = parent.parent().unwrap();
    let moved = f.root.join("relocated");
    fs::rename(old, &moved).unwrap();
    let parent = moved.join(parent.file_name().unwrap());
    let mut spec = f.spec("persisted.tsv");
    spec.source = EvidenceArtifactSource::Dataset {
        manifest: parent,
        query: DatasetQuery {
            selection: EvidenceSelection::WholeGenome,
            fields: fields(),
        },
        selection: ArtifactSelection::Stored,
        execution: EvidenceExecution::default(),
        limits: DatasetReadLimits::default(),
        artifacts: Vec::new(),
    };
    fs::remove_file(&f.bam).unwrap();
    fs::remove_file(&f.reference).unwrap();
    fs::remove_file(&f.bed).unwrap();
    let persisted = run_evidence_artifact(&mut Factory::default(), spec).unwrap();
    assert_eq!(
        fs::read(&native.output).unwrap(),
        fs::read(&persisted.output).unwrap()
    );
    assert_eq!(persisted.stats.record_visits, 0);
    let receipt = manifest(&persisted.manifest);
    assert_eq!(
        receipt.params["science.blake3"],
        manifest(&native.manifest).params["science.blake3"]
    );
    assert_eq!(receipt.inputs.len(), 8);
    assert!(receipt
        .inputs
        .iter()
        .all(|file| Path::new(&file.path).is_file()));
    assert!(receipt.params["command"].contains("--whole-dataset"));
    assert!(!receipt.params["command"].contains("--query-json"));
}

#[test]
fn factory_batch_finish_and_flush_failures_preserve_existing_artifact_and_receipt() {
    let _guard = SERIAL.lock().unwrap();
    let f = Fixture::new();
    for (index, failure) in [
        Failure::Create,
        Failure::Mismatch,
        Failure::Batch,
        Failure::Finish,
        Failure::Flush,
    ]
    .into_iter()
    .enumerate()
    {
        let mut spec = f.spec(&format!("failure{index}.tsv"));
        spec.output_policy = OutputPolicy::ReplaceAtomic;
        let receipt = f.root.join(format!("failure{index}.tsv.manifest.json"));
        fs::write(&spec.output, "old output").unwrap();
        fs::write(&receipt, "old receipt").unwrap();
        let output = spec.output.clone();
        assert!(run_evidence_artifact(
            &mut Factory {
                failure,
                unknown: false
            },
            spec
        )
        .is_err());
        assert_eq!(fs::read_to_string(output).unwrap(), "old output");
        assert_eq!(fs::read_to_string(receipt).unwrap(), "old receipt");
    }
    assert!(fs::read_dir(&f.root).unwrap().all(|entry| !entry
        .unwrap()
        .file_name()
        .to_string_lossy()
        .contains("rosalind-")));
}

#[test]
fn resource_breach_publishes_only_identified_partials() {
    let _guard = SERIAL.lock().unwrap();
    let f = Fixture::new();
    let spec = f.spec("limited.tsv");
    let output = spec.output.clone();
    let error = run_evidence_artifact(
        &mut Factory {
            failure: Failure::Resource,
            unknown: false,
        },
        spec,
    )
    .unwrap_err();
    assert_eq!(error.exit_code(), 4);
    let EvidenceArtifactError::Resource {
        partial: Some(partial),
        ..
    } = error
    else {
        panic!("missing resource partial")
    };
    assert!(!partial.completed);
    assert!(!output.exists());
    assert!(!f.root.join("limited.tsv.manifest.json").exists());
    assert!(partial.output.ends_with("limited.tsv.partial"));
    let receipt = manifest(&partial.manifest);
    assert_eq!(receipt.params["run_status"], "resource-failure");
    assert_eq!(
        receipt.params["artifact.output.0.role"],
        "partial-analyzer-output"
    );
    assert_eq!(
        receipt.outputs[0].blake3,
        blake3_file(&partial.output).unwrap()
    );
    assert_eq!(receipt.self_hash_ok(), Some(true));
    assert_eq!(receipt.measurement_hash_ok(), Some(true));
}

#[test]
fn cancellation_and_source_mutation_never_publish() {
    let _guard = SERIAL.lock().unwrap();
    for cancel_before in [true, false] {
        let f = Fixture::new();
        let token = CancellationToken::new();
        if cancel_before {
            token.cancel();
        }
        let mut spec = f.spec("cancel.tsv");
        spec.cancellation = Some(token.clone());
        let result = run_evidence_artifact(
            &mut Factory {
                failure: Failure::Cancel(token),
                unknown: false,
            },
            spec,
        );
        assert!(matches!(result, Err(EvidenceArtifactError::Cancelled)));
        assert!(!f.root.join("cancel.tsv").exists());
        assert!(!f.root.join("cancel.tsv.partial").exists());
    }
    let f = Fixture::new();
    let result = run_evidence_artifact(
        &mut Factory {
            failure: Failure::Mutate(f.bed.clone()),
            unknown: false,
        },
        f.spec("mutation.tsv"),
    );
    assert!(result.is_err());
    assert!(!f.root.join("mutation.tsv").exists());
    assert!(!f.root.join("mutation.tsv.partial").exists());
}

#[test]
fn tiny_budget_unknown_bound_invalid_claims_and_input_aliases_refuse_without_outputs() {
    let _guard = SERIAL.lock().unwrap();
    let f = Fixture::new();
    let mut spec = f.spec("refused.tsv");
    spec.enforcement = EnforcementMode::Cooperative;
    if let EvidenceArtifactSource::Native { request, .. } = &mut spec.source {
        request.execution.memory_budget_bytes = Some(1 << 20);
    }
    assert!(matches!(
        run_evidence_artifact(&mut Factory::default(), spec.clone()),
        Err(EvidenceArtifactError::Refused { .. })
    ));
    assert!(matches!(
        run_evidence_artifact(
            &mut Factory {
                unknown: true,
                ..Factory::default()
            },
            spec
        ),
        Err(EvidenceArtifactError::UnknownAnalyzerBound)
    ));
    let mut spec = f.spec("invalid.tsv");
    spec.invocation
        .options
        .insert("--fields".into(), "0".into());
    assert!(run_evidence_artifact(&mut Factory::default(), spec).is_err());
    let original = fs::read(&f.bam).unwrap();
    let mut spec = f.spec("ignored.tsv");
    spec.output = f.bam.clone();
    spec.output_policy = OutputPolicy::ReplaceAtomic;
    assert!(run_evidence_artifact(&mut Factory::default(), spec).is_err());
    assert_eq!(fs::read(&f.bam).unwrap(), original);
    assert!(!f.root.join("refused.tsv").exists());
    assert!(!f.root.join("invalid.tsv").exists());
}

#[test]
fn inline_queries_are_bounded_and_typed() {
    let _guard = SERIAL.lock().unwrap();
    let f = Fixture::new();
    let huge = " ".repeat(MAX_ARTIFACT_QUERY_BYTES + 1);
    assert!(parse_artifact_query(&huge)
        .unwrap_err()
        .to_string()
        .contains("file-backed"));
    assert!(parse_artifact_query(
        r#"{"version":1,"fields":128,"selection":{"intervals":[],"sites":null}}"#
    )
    .is_err());
    assert!(parse_artifact_query(
        r#"{"version":2,"fields":3,"selection":{"intervals":[],"sites":null}}"#
    )
    .is_err());
    assert!(parse_artifact_query(
        r#"{"version":1,"fields":3,"selection":{"intervals":[{"contig":0,"start":1,"end":2}],"sites":[{"contig":0,"position":2,"reference":65,"alternates":[67]}]}}"#
    ).is_err());
    assert!(parse_artifact_query(
        r#"{"version":1,"fields":3,"selection":{"intervals":[{"contig":0,"start":1,"end":2}],"sites":[{"contig":0,"position":1,"reference":65,"alternates":[67]}]}}"#
    ).is_ok());
    let mut spec = f.spec("large-query.tsv");
    if let EvidenceArtifactSource::Native { request, selection } = &mut spec.source {
        request.selection = EvidenceSelection::Sites(
            (0..1024)
                .map(|position| SnvSite {
                    contig: 0,
                    position,
                    reference: b'A',
                    alternates: vec![b'C'],
                })
                .collect(),
        );
        *selection = ArtifactSelection::Request;
    }
    assert!(run_evidence_artifact(&mut Factory::default(), spec)
        .unwrap_err()
        .to_string()
        .contains("file-backed"));
    assert!(!f.root.join("large-query.tsv").exists());
}

#[test]
fn large_file_selection_is_admitted_before_factory_construction() {
    let _guard = SERIAL.lock().unwrap();
    let f = Fixture::new();
    let variants = f.root.join("many.vcf");
    let mut text=String::from("##fileformat=VCFv4.2\n##contig=<ID=chr1,length=40000>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
    for position in 1..=25000 {
        text.push_str(&format!("chr1\t{position}\t.\tA\tC\t.\tPASS\t.\n"));
    }
    fs::write(&variants, text).unwrap();
    let mut spec = f.spec("large-file.tsv");
    spec.enforcement = EnforcementMode::Cooperative;
    if let EvidenceArtifactSource::Native { request, selection } = &mut spec.source {
        request.execution.memory_budget_bytes = Some(64 << 20);
        *selection =
            ArtifactSelection::Variants(variants, rosalind::variant_io::VariantLimits::default());
    }
    // A factory error would prove it ran too early. The size-based source
    // metadata reservation must refuse this budget before opening its sink.
    let result = run_evidence_artifact(
        &mut Factory {
            failure: Failure::Create,
            unknown: false,
        },
        spec,
    );
    assert!(
        matches!(result, Err(EvidenceArtifactError::Refused { .. })),
        "unexpected result: {result:?}"
    );
    assert!(!f.root.join("large-file.tsv").exists());
}
