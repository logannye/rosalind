//! Offline conformance for binaries using the exact evidence artifact runner.
use crate::conformance::ConformanceReport;
use crate::provenance::RunManifest;
use anyhow::{anyhow, Context, Result};
use rust_htslib::bam::{self, header::HeaderRecord, record::Cigar, record::CigarString};
use std::collections::BTreeMap;
use std::fs::{self, File};
use std::io::{Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::{Duration, Instant};

static NEXT: AtomicU64 = AtomicU64::new(0);
struct Workspace(PathBuf);
impl Drop for Workspace {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

/// Exercise standard `run` operands, canonical budget/tile output, managed receipt
/// identity, native/persisted equality, resource partials, mutation rejection and
/// explicit relocated replay. Reducer-specific scientific oracles remain the
/// analyzer author's tests; this harness does not assume a particular TSV schema.
pub fn conform_evidence_analyzer(binary: &Path) -> Result<ConformanceReport> {
    let binary = fs::canonicalize(binary).context("cannot locate external evidence analyzer")?;
    let rosalind = std::env::current_exe()?;
    let root = Workspace(std::env::temp_dir().join(format!(
        "rosalind-evidence-conformance-{}-{}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    )));
    fs::create_dir(&root.0)?;
    let fixture = fixture(&root.0)?;
    let mut report = ConformanceReport {
        passed: false,
        checks: BTreeMap::new(),
        failures: Vec::new(),
    };
    let mut artifacts = Vec::new();
    for (budget, tile) in [(128, 1), (256, 137), (512, 16384)] {
        let output = root.0.join(format!("native-{budget}.tsv"));
        let mut command = native(&binary, &fixture, &output);
        command.args([
            "--memory-budget-mb",
            &budget.to_string(),
            "--tile-bases",
            &tile.to_string(),
            "--enforce",
        ]);
        let result = run(&mut command, &root.0)?;
        check(
            &mut report,
            &format!("admitted_{budget}MiB_tile_{tile}"),
            result.code == Some(0) && output.is_file(),
            result.detail(),
        );
        artifacts.push(output);
    }
    let first = &artifacts[0];
    let equal = crate::provenance::blake3_file(first)
        .ok()
        .is_some_and(|bytes| {
            artifacts
                .iter()
                .skip(1)
                .all(|p| crate::provenance::blake3_file(p).is_ok_and(|next| next == bytes))
        });
    check(
        &mut report,
        "canonical_budget_tile_byte_equality",
        equal,
        "outputs differ or are missing".into(),
    );
    let first_receipt = receipt(&sidecar(first, ".manifest.json"));
    check(
        &mut report,
        "external_receipt_identity_and_integrity",
        first_receipt.as_ref().is_ok_and(|m| {
            m.self_hash_ok() == Some(true)
                && m.measurement_hash_ok() == Some(true)
                && m.params
                    .get("replay.kind")
                    .is_some_and(|v| v == "external-analyzer")
                && m.params.get("run_status").is_some_and(|v| v == "completed")
                && m.params.contains_key("evidence.artifact_semantics")
        }),
        first_receipt
            .as_ref()
            .err()
            .map(ToString::to_string)
            .unwrap_or_else(|| "incomplete artifact identity".into()),
    );
    if let Ok(first_receipt) = &first_receipt {
        let same_science = artifacts.iter().skip(1).all(|p| {
            receipt(&sidecar(p, ".manifest.json")).is_ok_and(|m| {
                first_receipt.params.contains_key("science.blake3")
                    && m.params.get("science.blake3") == first_receipt.params.get("science.blake3")
            })
        });
        check(
            &mut report,
            "execution_settings_do_not_change_scientific_identity",
            same_science,
            "scientific identity depends on budget or tile".into(),
        );
    }
    let low = root.0.join("refused.tsv");
    let result = run(
        native(&binary, &fixture, &low).args(["--memory-budget-mb", "1", "--enforce"]),
        &root.0,
    )?;
    check(
        &mut report,
        "insufficient_budget_has_no_success_artifact",
        matches!(result.code, Some(3 | 4)) && !low.exists(),
        result.detail(),
    );
    let limited = root.0.join("limited.tsv");
    let result = run(
        native(&binary, &fixture, &limited).args([
            "--max-read-len",
            "1",
            "--memory-budget-mb",
            "128",
            "--enforce",
        ]),
        &root.0,
    )?;
    let partial_receipt = receipt(&sidecar(&limited, ".manifest.json.partial"));
    check(
        &mut report,
        "declared_record_limit_has_identified_partial",
        result.code == Some(4)
            && !limited.exists()
            && sidecar(&limited, ".partial").is_file()
            && partial_receipt.is_ok_and(|m| {
                m.self_hash_ok() == Some(true)
                    && m.params.get("run_status").is_some_and(|v| v != "completed")
            }),
        result.detail(),
    );
    let first_bytes = crate::provenance::blake3_file(first).ok();
    let result = run(&mut native(&binary, &fixture, first), &root.0)?;
    check(
        &mut report,
        "safe_output_collision",
        result.code == Some(2) && crate::provenance::blake3_file(first).ok() == first_bytes,
        result.detail(),
    );
    let result = run(native(&binary, &fixture, first).arg("--force"), &root.0)?;
    check(
        &mut report,
        "explicit_atomic_replacement",
        result.code == Some(0) && crate::provenance::blake3_file(first).ok() == first_bytes,
        result.detail(),
    );

    let changed = root.0.join("changed-profile.tsv");
    let result = run(
        native(&binary, &fixture, &changed).args(["--mapq-threshold", "61"]),
        &root.0,
    )?;
    let changed_receipt = receipt(&sidecar(&changed, ".manifest.json"));
    check(
        &mut report,
        "scientific_filter_changes_are_recorded",
        result.code == Some(0)
            && first_receipt
                .as_ref()
                .ok()
                .zip(changed_receipt.as_ref().ok())
                .is_some_and(|(a, b)| {
                    a.params.contains_key("science.blake3")
                        && a.params.get("science.blake3") != b.params.get("science.blake3")
                }),
        result.detail(),
    );

    let relocated = root.0.join("relocated");
    fs::create_dir(&relocated)?;
    if let Ok(manifest) = &first_receipt {
        for (index, input) in manifest.inputs.iter().enumerate() {
            fs::copy(&input.path, relocated.join(format!("input-{index}")))?;
        }
        let _hidden = HiddenInputs::hide(manifest, &root.0)?;
        let result = run(
            replay(
                &rosalind,
                &sidecar(first, ".manifest.json"),
                &relocated,
                None,
            )
            .arg("--dry-run"),
            &root.0,
        )?;
        check(
            &mut report,
            "external_replay_requires_explicit_binary",
            result.code != Some(0),
            result.detail(),
        );
        let result = run(
            &mut replay(
                &rosalind,
                &sidecar(first, ".manifest.json"),
                &relocated,
                Some(&binary),
            ),
            &root.0,
        )?;
        check(
            &mut report,
            "relocated_native_byte_replay",
            result.code == Some(0) && result.stdout.contains("REPRODUCED"),
            result.detail(),
        );
    }

    let custom = root.0.join("encoded.custom");
    let result = run(
        native(&binary, &fixture, &custom).args([
            "--memory-budget-bytes",
            "134217735",
            "--enforce",
        ]),
        &root.0,
    )?;
    let custom_receipt = receipt(&sidecar(&custom, ".manifest.json"));
    check(
        &mut report,
        "exact_byte_budget_and_custom_output",
        result.code == Some(0)
            && crate::provenance::blake3_file(&custom).ok() == first_bytes
            && custom_receipt.as_ref().is_ok_and(|m| {
                let checked = crate::provenance::verify_receipt(
                    &m.to_canonical_json(),
                    &crate::provenance::VerifyOpts::default(),
                );
                checked.ok && m.memory_budget_bytes().ok().flatten() == Some(134217735)
            }),
        result.detail(),
    );
    let result = run(
        &mut replay(
            &rosalind,
            &sidecar(&custom, ".manifest.json"),
            &relocated,
            Some(&binary),
        ),
        &root.0,
    )?;
    check(
        &mut report,
        "custom_output_has_physical_byte_replay",
        result.code == Some(0) && result.stdout.contains("REPRODUCED"),
        result.detail(),
    );

    // Persist full capabilities, then use the external factory's projection.
    let raw = root.0.join("raw.arrow");
    let cache = root.0.join("cache");
    let mut create = Command::new(&rosalind);
    create
        .args(["analyze", "evidence", "--reference"])
        .arg(&fixture.reference)
        .arg("--alignments")
        .arg(&fixture.alignment)
        .arg("--regions")
        .arg(&fixture.regions)
        .args([
            "--fields",
            "all-supported",
            "--format",
            "arrow-ipc",
            "--cache-dir",
        ])
        .arg(&cache)
        .arg("-o")
        .arg(&raw);
    let created = run(&mut create, &root.0)?;
    check(
        &mut report,
        "verified_dataset_fixture",
        created.code == Some(0),
        created.detail(),
    );
    if created.code == Some(0) {
        let raw_receipt = receipt(&sidecar(&raw, ".manifest.json"))?;
        let portable = PathBuf::from(
            raw_receipt
                .measurements
                .get("execution.evidence_dataset_manifest")
                .ok_or_else(|| anyhow!("dataset fixture did not publish portable metadata"))?,
        );
        let moved = relocated.join("portable");
        fs::rename(portable.parent().unwrap(), &moved)?;
        let manifest = moved.join("evidence-dataset.manifest.json");
        // Physical source deletion proves the external dataset path is independent.
        fs::remove_file(&fixture.alignment)?;
        fs::remove_file(&fixture.reference)?;
        let derived = root.0.join("persisted.tsv");
        let mut query = Command::new(&binary);
        query
            .args(["run", "--dataset"])
            .arg(&manifest)
            .args(["--whole-dataset", "--output"])
            .arg(&derived);
        let result = run(&mut query, &root.0)?;
        check(
            &mut report,
            "persisted_source_matches_native_without_alignments",
            result.code == Some(0) && crate::provenance::blake3_file(&derived).ok() == first_bytes,
            result.detail(),
        );
        let result = run(
            &mut replay(
                &rosalind,
                &sidecar(&derived, ".manifest.json"),
                &relocated,
                Some(&binary),
            ),
            &root.0,
        )?;
        check(
            &mut report,
            "relocated_persisted_byte_replay",
            result.code == Some(0) && result.stdout.contains("REPRODUCED"),
            result.detail(),
        );
        let descriptor: serde_json::Value =
            serde_json::from_slice(&read_metadata(&moved.join("dataset.descriptor.json"))?)?;
        if let Some(part) = descriptor["partitions"]
            .as_array()
            .and_then(|parts| parts.first())
        {
            let path = moved.join(
                part["arrow"]["path"]
                    .as_str()
                    .ok_or_else(|| anyhow!("fixture Arrow path absent"))?,
            );
            let mut file = fs::OpenOptions::new().read(true).write(true).open(path)?;
            let offset = file.metadata()?.len() / 2;
            file.seek(SeekFrom::Start(offset))?;
            let mut byte = [0u8];
            file.read_exact(&mut byte)?;
            file.seek(SeekFrom::Start(offset))?;
            file.write_all(&[byte[0] ^ 1])?;
            file.sync_all()?;
            let old_receipt = crate::provenance::blake3_file(&sidecar(&derived, ".manifest.json"))?;
            let result = run(query.arg("--force"), &root.0)?;
            check(
                &mut report,
                "corrupt_source_cannot_replace_existing_result",
                result.code == Some(5)
                    && crate::provenance::blake3_file(&derived).ok() == first_bytes
                    && crate::provenance::blake3_file(&sidecar(&derived, ".manifest.json"))?
                        == old_receipt,
                result.detail(),
            );
        }
    }
    report.passed = report.checks.values().all(|passed| *passed);
    Ok(report)
}

struct HiddenInputs(Vec<(PathBuf, PathBuf)>);
impl HiddenInputs {
    fn hide(manifest: &RunManifest, root: &Path) -> Result<Self> {
        let mut hidden = Self(Vec::new());
        let paths = manifest
            .inputs
            .iter()
            .map(|input| PathBuf::from(&input.path))
            .collect::<std::collections::BTreeSet<_>>();
        for (index, source) in paths.into_iter().enumerate() {
            let destination = root.join(format!("temporarily-hidden-input-{index}"));
            fs::rename(&source, &destination)?;
            hidden.0.push((source, destination));
        }
        Ok(hidden)
    }
}
impl Drop for HiddenInputs {
    fn drop(&mut self) {
        for (original, hidden) in self.0.iter().rev() {
            let _ = fs::rename(hidden, original);
        }
    }
}

struct Fixture {
    reference: PathBuf,
    alignment: PathBuf,
    regions: PathBuf,
}
fn fixture(root: &Path) -> Result<Fixture> {
    let reference = root.join("reference.fa");
    fs::write(&reference, format!(">chr1\n{}\n", "A".repeat(32800)))?;
    fs::write(
        root.join("reference.fa.fai"),
        "chr1\t32800\t6\t32800\t32801\n",
    )?;
    let alignment = root.join("reads.bam");
    let mut header = bam::Header::new();
    header.push_record(
        HeaderRecord::new(b"HD")
            .push_tag(b"VN", "1.6")
            .push_tag(b"SO", "coordinate"),
    );
    header.push_record(
        HeaderRecord::new(b"SQ")
            .push_tag(b"SN", "chr1")
            .push_tag(b"LN", 32800),
    );
    let mut writer = bam::Writer::from_path(&alignment, &header, bam::Format::Bam)?;
    for start in [0, 1020, 16370, 32750] {
        for index in 0..12 {
            let mut read = bam::Record::new();
            read.set(
                format!("r{start}-{index}").as_bytes(),
                Some(&CigarString(vec![Cigar::Match(30)])),
                &[if index < 8 { b'A' } else { b'C' }; 30],
                &[30; 30],
            );
            read.set_tid(0);
            read.set_pos(start);
            read.set_mapq(60);
            read.set_flags(if index % 2 == 0 { 0 } else { 16 });
            writer.write(&read)?;
        }
    }
    drop(writer);
    bam::index::build(&alignment, None::<&PathBuf>, bam::index::Type::Bai, 1)?;
    let regions = root.join("selection.bed");
    fs::write(
        &regions,
        "chr1\t0\t2050\nchr1\t16370\t16410\nchr1\t32750\t32800\n",
    )?;
    Ok(Fixture {
        reference,
        alignment,
        regions,
    })
}
fn native(binary: &Path, fixture: &Fixture, output: &Path) -> Command {
    let mut command = Command::new(binary);
    command
        .args(["run", "--reference"])
        .arg(&fixture.reference)
        .arg("--alignments")
        .arg(&fixture.alignment)
        .arg("--regions")
        .arg(&fixture.regions)
        .arg("--output")
        .arg(output);
    command
}
fn replay(rosalind: &Path, manifest: &Path, inputs: &Path, binary: Option<&Path>) -> Command {
    let mut command = Command::new(rosalind);
    command
        .args(["reproduce", "--manifest"])
        .arg(manifest)
        .arg("--inputs")
        .arg(inputs)
        .args(["--json", "--no-attest"]);
    if let Some(binary) = binary {
        command.arg("--binary").arg(binary);
    }
    command
}
fn sidecar(path: &Path, suffix: &str) -> PathBuf {
    let mut name = path.as_os_str().to_os_string();
    name.push(suffix);
    PathBuf::from(name)
}
fn receipt(path: &Path) -> Result<RunManifest> {
    let bytes = read_metadata(path)?;
    RunManifest::from_canonical_json(std::str::from_utf8(&bytes)?).map_err(|e| anyhow!(e))
}
fn read_metadata(path: &Path) -> Result<Vec<u8>> {
    const LIMIT: u64 = 32 << 20;
    let file = File::open(path)?;
    if file.metadata()?.len() > LIMIT {
        anyhow::bail!("metadata exceeds32MiB: {}", path.display());
    }
    let mut bytes = Vec::new();
    file.take(LIMIT + 1).read_to_end(&mut bytes)?;
    if bytes.len() as u64 > LIMIT {
        anyhow::bail!("metadata grew beyond32MiB: {}", path.display());
    }
    Ok(bytes)
}
fn check(report: &mut ConformanceReport, name: &str, passed: bool, detail: String) {
    report.checks.insert(name.into(), passed);
    if !passed {
        report.failures.push(format!("{name}: {detail}"));
    }
}
struct ChildResult {
    code: Option<i32>,
    stdout: String,
    stderr: String,
    timed_out: bool,
}
impl ChildResult {
    fn detail(&self) -> String {
        format!(
            "exit {:?}, timed_out {}: {} {}",
            self.code, self.timed_out, self.stderr, self.stdout
        )
    }
}
fn run(command: &mut Command, root: &Path) -> Result<ChildResult> {
    let id = NEXT.fetch_add(1, Ordering::Relaxed);
    let stdout = root.join(format!("stdout-{id}.log"));
    let stderr = root.join(format!("stderr-{id}.log"));
    command
        .stdout(Stdio::from(File::create(&stdout)?))
        .stderr(Stdio::from(File::create(&stderr)?));
    #[cfg(unix)]
    {
        use std::os::unix::process::CommandExt;
        command.process_group(0);
    }
    let mut child = command.spawn()?;
    let start = Instant::now();
    let mut timed_out = false;
    let status = loop {
        if let Some(status) = child.try_wait()? {
            break status;
        }
        if start.elapsed() > Duration::from_secs(60) {
            timed_out = true;
            #[cfg(unix)]
            {
                // SAFETY: process_group(0) gives only this owned child and its
                // descendants a new group; no unrelated processes share its id.
                unsafe {
                    libc::kill(-(child.id() as i32), libc::SIGKILL);
                }
            }
            let _ = child.kill();
            break child.wait()?;
        }
        std::thread::sleep(Duration::from_millis(10));
    };
    Ok(ChildResult {
        code: status.code(),
        stdout: tail(&stdout)?,
        stderr: tail(&stderr)?,
        timed_out,
    })
}
fn tail(path: &Path) -> Result<String> {
    let mut f = File::open(path)?;
    let n = f.metadata()?.len();
    f.seek(SeekFrom::Start(n.saturating_sub(32768)))?;
    let mut bytes = Vec::new();
    f.take(32768).read_to_end(&mut bytes)?;
    Ok(String::from_utf8_lossy(&bytes).into_owned())
}
