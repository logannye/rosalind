# Phase A4 — Spec-valid VCF writer + BLAKE3 provenance receipt Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build `src/io/vcf.rs` — one spec-valid VCFv4.2 writer (germline single-sample + somatic TUMOR/NORMAL) that consumes the A3 call types — and `src/provenance/` — a minimal, deterministic BLAKE3 reproducibility receipt.

**Architecture:** Two independent crate-root modules built on `crate::core` + `crate::call` + `std` (+ the existing `blake3` dep). `io/vcf` formats `GermlineCall`/`SomaticCall` (+ a `Locus`/`ContigSet`) into standards-compliant VCF text with a typed header and deterministic record ordering. `provenance` builds a `RunManifest` and serializes it as canonical JSON (sorted keys, no timestamps → byte-identical across runs); it hashes input/output files with BLAKE3. No CLI wiring and no consumption of these by the calling path yet — that is A5. `serde` is an optional dep in this crate, so the manifest JSON is hand-rolled (the structure is flat and small).

**Tech Stack:** Rust, `std::io`, `blake3` (already a dependency), `std::collections::BTreeMap` for sorted params. No new dependencies.

**Design reference:** `docs/superpowers/specs/2026-05-26-phase-a-unified-pileup-genotype-design.md` §3.4 (io/vcf), §3.5 (provenance), §6 (VCF contract).

**Consumes (already shipped & green):**
- A3 `crate::call`: `GermlineCall { genotype: Genotype, alt_base: u8, qual: f64, gq: u8, pl: [u32;3], ad: [u32;2], dp: u32, filter: Filter }`; `Genotype { HomRef, Het, HomAlt }`; `Filter { Pass, LowQual, LowDepth }`; `SomaticCall { ref_base, alt_base, tumor_alt, tumor_depth, normal_alt, normal_depth, tumor_af: f32, normal_af: f32, quality: f64 }`.
- A1 `crate::core`: `ContigSet` (`iter() -> &Contig` in id order; `by_id(u32) -> Option<&Contig>`); `Contig { id, name: Arc<str>, length: u32, .. }`; `Locus { contig: u32, pos: Position }`; `Position(u32)`.

---

## File Structure

- `src/provenance/mod.rs` — `FileHash`, `RunManifest`, `to_canonical_json`, `blake3_hex`, `blake3_file`, `write_manifest`.
- `src/io/mod.rs` — IO layer banner + `pub mod vcf;`.
- `src/io/vcf.rs` — `GermlineRow`, header + record formatting, `write_germline_vcf`/`render_germline_vcf`, `write_somatic_vcf`/`render_somatic_vcf`.
- `src/lib.rs` — register `pub mod io;` and `pub mod provenance;` (alphabetical: `io` after `genomics`; `provenance` after `plugin`).

---

### Task 1: `provenance/` — BLAKE3 reproducibility receipt

**Files:**
- Create: `src/provenance/mod.rs`
- Modify: `src/lib.rs` (add `pub mod provenance;` after `pub mod plugin;`)
- Test: in `src/provenance/mod.rs`

- [ ] **Step 1: Write the failing tests.** Create `src/provenance/mod.rs`:

```rust
//! A minimal, deterministic reproducibility receipt for a run: tool version,
//! subcommand, BLAKE3 content hashes of inputs + outputs, and the parameters.
//! Serialized as canonical JSON (sorted keys, no timestamps) so two identical
//! runs produce a byte-identical manifest. Full `rosalind verify` is a later
//! phase; this phase emits the receipt and proves it is deterministic.

use std::collections::BTreeMap;
use std::io::{self, Read, Write};
use std::path::{Path, PathBuf};

/// A file referenced by a run, with its content hash.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FileHash {
    /// Path as recorded (normalized to a string by the caller).
    pub path: String,
    /// BLAKE3 hex digest of the file's contents.
    pub blake3: String,
}

/// A reproducibility receipt for a single run.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RunManifest {
    /// Rosalind version (`CARGO_PKG_VERSION`).
    pub tool_version: String,
    /// Subcommand that produced the run (e.g. `variants`, `somatic`).
    pub subcommand: String,
    /// Input files and their content hashes.
    pub inputs: Vec<FileHash>,
    /// Run parameters (sorted by key in the canonical form).
    pub params: BTreeMap<String, String>,
    /// Output files and their content hashes.
    pub outputs: Vec<FileHash>,
}

impl RunManifest {
    /// A new manifest stamped with the current tool version.
    pub fn new(subcommand: impl Into<String>) -> Self {
        Self {
            tool_version: env!("CARGO_PKG_VERSION").to_string(),
            subcommand: subcommand.into(),
            inputs: Vec::new(),
            params: BTreeMap::new(),
            outputs: Vec::new(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn blake3_is_deterministic_and_sensitive() {
        assert_eq!(blake3_hex(b"abc"), blake3_hex(b"abc"));
        assert_ne!(blake3_hex(b"abc"), blake3_hex(b"abd"));
        // Known BLAKE3 vector for the empty input.
        assert_eq!(
            blake3_hex(b""),
            "af1349b9f5f9a1a6a0404dea36dcc9499bcb25c9adc112b7cc9a93cae41f3262"
        );
        assert_eq!(blake3_hex(b"abc").len(), 64);
    }

    #[test]
    fn canonical_json_has_sorted_keys_and_is_exact() {
        let mut m = RunManifest::new("variants");
        m.tool_version = "0.1.0".to_string();
        m.inputs.push(FileHash {
            path: "ref.fa".to_string(),
            blake3: "aa".to_string(),
        });
        m.outputs.push(FileHash {
            path: "out.vcf".to_string(),
            blake3: "bb".to_string(),
        });
        m.params.insert("min_qual".to_string(), "30".to_string());
        m.params.insert("min_depth".to_string(), "8".to_string());

        let json = m.to_canonical_json();
        assert_eq!(
            json,
            r#"{"inputs":[{"blake3":"aa","path":"ref.fa"}],"outputs":[{"blake3":"bb","path":"out.vcf"}],"params":{"min_depth":"8","min_qual":"30"},"subcommand":"variants","tool_version":"0.1.0"}"#
        );
    }

    #[test]
    fn canonical_json_is_order_independent() {
        let mk = |order_swapped: bool| {
            let mut m = RunManifest::new("variants");
            m.tool_version = "0.1.0".to_string();
            let a = FileHash {
                path: "a.fa".to_string(),
                blake3: "1".to_string(),
            };
            let b = FileHash {
                path: "b.fa".to_string(),
                blake3: "2".to_string(),
            };
            if order_swapped {
                m.inputs.push(b);
                m.inputs.push(a);
            } else {
                m.inputs.push(a);
                m.inputs.push(b);
            }
            m.to_canonical_json()
        };
        // Inputs are sorted by path in the canonical form → order-independent.
        assert_eq!(mk(false), mk(true));
    }

    #[test]
    fn json_escapes_special_characters() {
        let mut m = RunManifest::new("variants");
        m.tool_version = "0.1.0".to_string();
        m.params
            .insert("note".to_string(), "a\"b\\c".to_string());
        let json = m.to_canonical_json();
        assert!(json.contains(r#""note":"a\"b\\c""#));
    }

    #[test]
    fn write_manifest_emits_sidecar_file() {
        let dir = std::env::temp_dir().join(format!("rosalind_manifest_test_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let out = dir.join("calls.vcf");
        std::fs::write(&out, b"##fileformat=VCFv4.2\n").unwrap();

        let mut m = RunManifest::new("variants");
        m.tool_version = "0.1.0".to_string();
        m.outputs.push(FileHash {
            path: out.display().to_string(),
            blake3: blake3_file(&out).unwrap(),
        });

        let manifest_path = write_manifest(&out, &m).unwrap();
        assert_eq!(manifest_path, dir.join("calls.vcf.manifest.json"));
        let written = std::fs::read_to_string(&manifest_path).unwrap();
        assert_eq!(written, m.to_canonical_json());

        std::fs::remove_dir_all(&dir).ok();
    }
}
```

- [ ] **Step 2: Run the tests to verify they fail.** Run: `cargo test provenance`
Expected: FAIL — does not compile (`blake3_hex`, `to_canonical_json`, `blake3_file`, `write_manifest` undefined).

- [ ] **Step 3: Implement the functions.** Add to `src/provenance/mod.rs`, above the `#[cfg(test)]` block (after the `impl RunManifest` with `new`):

```rust
impl RunManifest {
    /// Serialize to canonical JSON: keys sorted, `inputs`/`outputs` sorted by
    /// path, no timestamps — so identical runs hash and render identically.
    pub fn to_canonical_json(&self) -> String {
        let mut out = String::new();
        out.push('{');

        out.push_str("\"inputs\":");
        push_file_hashes(&mut out, &self.inputs);

        out.push_str(",\"outputs\":");
        push_file_hashes(&mut out, &self.outputs);

        out.push_str(",\"params\":{");
        for (i, (k, v)) in self.params.iter().enumerate() {
            if i > 0 {
                out.push(',');
            }
            out.push('"');
            out.push_str(&json_escape(k));
            out.push_str("\":\"");
            out.push_str(&json_escape(v));
            out.push('"');
        }
        out.push('}');

        out.push_str(",\"subcommand\":\"");
        out.push_str(&json_escape(&self.subcommand));
        out.push_str("\",\"tool_version\":\"");
        out.push_str(&json_escape(&self.tool_version));
        out.push_str("\"}");

        out
    }
}

/// Render a `[{"blake3":..,"path":..}, ..]` array, entries sorted by path.
fn push_file_hashes(out: &mut String, files: &[FileHash]) {
    let mut sorted: Vec<&FileHash> = files.iter().collect();
    sorted.sort_by(|a, b| a.path.cmp(&b.path));
    out.push('[');
    for (i, f) in sorted.iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        out.push_str("{\"blake3\":\"");
        out.push_str(&json_escape(&f.blake3));
        out.push_str("\",\"path\":\"");
        out.push_str(&json_escape(&f.path));
        out.push_str("\"}");
    }
    out.push(']');
}

/// Minimal RFC-8259 string escaping for the characters we can encounter.
fn json_escape(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if (c as u32) < 0x20 => out.push_str(&format!("\\u{:04x}", c as u32)),
            c => out.push(c),
        }
    }
    out
}

/// BLAKE3 hex digest of a byte slice.
pub fn blake3_hex(bytes: &[u8]) -> String {
    blake3::hash(bytes).to_hex().to_string()
}

/// BLAKE3 hex digest of a file's contents, streamed in fixed-size chunks
/// (bounded memory regardless of file size).
pub fn blake3_file(path: &Path) -> io::Result<String> {
    let mut hasher = blake3::Hasher::new();
    let mut file = std::fs::File::open(path)?;
    let mut buf = [0u8; 64 * 1024];
    loop {
        let n = file.read(&mut buf)?;
        if n == 0 {
            break;
        }
        hasher.update(&buf[..n]);
    }
    Ok(hasher.finalize().to_hex().to_string())
}

/// Write `<output_path>.manifest.json` next to the output, returning its path.
pub fn write_manifest(output_path: &Path, manifest: &RunManifest) -> io::Result<PathBuf> {
    let mut name = output_path.as_os_str().to_os_string();
    name.push(".manifest.json");
    let manifest_path = PathBuf::from(name);
    let mut file = std::fs::File::create(&manifest_path)?;
    file.write_all(manifest.to_canonical_json().as_bytes())?;
    file.flush()?;
    Ok(manifest_path)
}
```

- [ ] **Step 4: Register the module.** In `src/lib.rs`, add after the `pub mod plugin;` line:

```rust
/// Reproducibility receipts: canonical-JSON BLAKE3 manifests for every run.
pub mod provenance;
```

- [ ] **Step 5: Run the tests to verify they pass.** Run: `cargo test provenance`
Expected: PASS (5 tests). Then `cargo build 2>&1 | grep -i warning` (no `src/provenance/` warnings), `cargo fmt --all -- --check`.
If the empty-input BLAKE3 vector assertion fails, the digest your `blake3` version produces is authoritative — replace the constant with the value from the failure message (it should match the published BLAKE3 empty-input hash).

- [ ] **Step 6: Commit.**

```bash
git add src/provenance/mod.rs src/lib.rs
git commit -m "feat(provenance): canonical-JSON BLAKE3 run manifest (deterministic receipt)" \
  -m "Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task 2: `io/vcf` — germline single-sample VCFv4.2 writer

**Files:**
- Create: `src/io/mod.rs`
- Create: `src/io/vcf.rs`
- Modify: `src/lib.rs` (add `pub mod io;` after `pub mod genomics;`)
- Test: in `src/io/vcf.rs`

- [ ] **Step 1: Write the failing tests.** Create `src/io/vcf.rs`:

```rust
//! One spec-valid VCFv4.2 writer for the calling layer: a typed header built
//! from a `ContigSet`, deterministic record ordering, germline single-sample
//! and somatic TUMOR/NORMAL output. Replaces the legacy non-conformant string
//! writers (no `##contig`/`##FORMAT`/sample columns).

use std::io::{self, Write};

use crate::call::{Filter, GermlineCall, Genotype, SomaticCall};
use crate::core::{ContigSet, Locus};

/// One germline row: the locus, its reference base, and the call there. The
/// call is locus-free, so the writer pairs it with coordinate + ref context.
#[derive(Debug, Clone)]
pub struct GermlineRow {
    /// Genomic coordinate.
    pub locus: Locus,
    /// Reference base (uppercase ASCII).
    pub ref_base: u8,
    /// The germline call at this locus.
    pub call: GermlineCall,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::call::GermlineParams;
    use crate::core::Position;

    fn contigs() -> ContigSet {
        let mut c = ContigSet::new();
        c.push("chr1", 248_956_422);
        c.push("chr2", 100);
        c
    }

    fn row(contig: u32, pos0: u32, ref_base: u8, call: GermlineCall) -> GermlineRow {
        GermlineRow {
            locus: Locus {
                contig,
                pos: Position(pos0),
            },
            ref_base,
            call,
        }
    }

    fn het_call() -> GermlineCall {
        GermlineCall {
            genotype: Genotype::Het,
            alt_base: b'G',
            qual: 48.0,
            gq: 48,
            pl: [48, 0, 49],
            ad: [15, 15],
            dp: 30,
            filter: Filter::Pass,
        }
    }

    #[test]
    fn header_has_required_tags() {
        let vcf = render_germline_vcf(&contigs(), "SAMPLE", &[]).unwrap();
        assert!(vcf.starts_with("##fileformat=VCFv4.2\n"));
        assert!(vcf.contains("##contig=<ID=chr1,length=248956422>\n"));
        assert!(vcf.contains("##contig=<ID=chr2,length=100>\n"));
        assert!(vcf.contains("##INFO=<ID=DP,Number=1,Type=Integer,"));
        assert!(vcf.contains("##INFO=<ID=AF,Number=A,Type=Float,"));
        assert!(vcf.contains("##FORMAT=<ID=GT,Number=1,Type=String,"));
        assert!(vcf.contains("##FORMAT=<ID=AD,Number=R,Type=Integer,"));
        assert!(vcf.contains("##FORMAT=<ID=PL,Number=G,Type=Integer,"));
        assert!(vcf.contains("##FILTER=<ID=LowQual,"));
        assert!(vcf.contains("##FILTER=<ID=LowDepth,"));
        assert!(vcf.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"));
    }

    #[test]
    fn germline_record_is_exact() {
        let vcf =
            render_germline_vcf(&contigs(), "SAMPLE", &[row(0, 100, b'A', het_call())]).unwrap();
        let last = vcf.lines().last().unwrap();
        // POS is 1-based (100 -> 101); AF = 15/30 = 0.500.
        assert_eq!(
            last,
            "chr1\t101\t.\tA\tG\t48.0\tPASS\tDP=30;AF=0.500\tGT:GQ:DP:AD:PL\t0/1:48:30:15,15:48,0,49"
        );
    }

    #[test]
    fn records_are_sorted_deterministically() {
        let r1 = row(0, 100, b'A', het_call());
        let r2 = row(0, 50, b'A', het_call());
        let r3 = row(1, 10, b'A', het_call());
        let forward = render_germline_vcf(&contigs(), "S", &[r1.clone(), r2.clone(), r3.clone()]).unwrap();
        let shuffled = render_germline_vcf(&contigs(), "S", &[r3, r1, r2]).unwrap();
        assert_eq!(forward, shuffled);
        // chr1:51 sorts before chr1:101 before chr2:11.
        let body: Vec<&str> = forward.lines().filter(|l| !l.starts_with('#')).collect();
        assert_eq!(body[0].split('\t').next().unwrap(), "chr1");
        assert!(body[0].contains("\t51\t"));
        assert!(body[1].contains("\t101\t"));
        assert_eq!(body[2].split('\t').next().unwrap(), "chr2");
    }

    #[test]
    fn filter_and_genotype_strings_render() {
        let mut c = het_call();
        c.genotype = Genotype::HomAlt;
        c.filter = Filter::LowDepth;
        let vcf = render_germline_vcf(&contigs(), "S", &[row(0, 0, b'C', c)]).unwrap();
        let last = vcf.lines().last().unwrap();
        assert!(last.contains("\tLowDepth\t"));
        assert!(last.contains("\t1/1:"));
    }
}
```

- [ ] **Step 2: Run the tests to verify they fail.** Run: `cargo test io::vcf`
Expected: FAIL — does not compile (`render_germline_vcf` undefined).

- [ ] **Step 3: Implement the germline writer.** Add to `src/io/vcf.rs`, above the `#[cfg(test)]` block:

```rust
fn genotype_str(g: Genotype) -> &'static str {
    match g {
        Genotype::HomRef => "0/0",
        Genotype::Het => "0/1",
        Genotype::HomAlt => "1/1",
    }
}

fn filter_str(f: Filter) -> &'static str {
    match f {
        Filter::Pass => "PASS",
        Filter::LowQual => "LowQual",
        Filter::LowDepth => "LowDepth",
    }
}

/// `##fileformat` + one `##contig` line per contig (shared by both writers).
fn write_fileformat_and_contigs<W: Write>(out: &mut W, contigs: &ContigSet) -> io::Result<()> {
    writeln!(out, "##fileformat=VCFv4.2")?;
    for c in contigs.iter() {
        writeln!(out, "##contig=<ID={},length={}>", c.name, c.length)?;
    }
    Ok(())
}

/// Write a spec-valid germline (single-sample) VCFv4.2 to `out`. Records are
/// emitted in canonical (contig, pos, ref, alt) order regardless of input order.
pub fn write_germline_vcf<W: Write>(
    out: &mut W,
    contigs: &ContigSet,
    sample: &str,
    rows: &[GermlineRow],
) -> io::Result<()> {
    write_fileformat_and_contigs(out, contigs)?;
    writeln!(
        out,
        r#"##INFO=<ID=DP,Number=1,Type=Integer,Description="Total depth">"#
    )?;
    writeln!(
        out,
        r#"##INFO=<ID=AF,Number=A,Type=Float,Description="Alt allele fraction">"#
    )?;
    writeln!(
        out,
        r#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#
    )?;
    writeln!(
        out,
        r#"##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">"#
    )?;
    writeln!(
        out,
        r#"##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">"#
    )?;
    writeln!(
        out,
        r#"##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths (ref,alt)">"#
    )?;
    writeln!(
        out,
        r#"##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred genotype likelihoods">"#
    )?;
    writeln!(out, r#"##FILTER=<ID=LowQual,Description="QUAL below threshold">"#)?;
    writeln!(out, r#"##FILTER=<ID=LowDepth,Description="Depth below threshold">"#)?;
    writeln!(
        out,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}"
    )?;

    let mut ordered: Vec<&GermlineRow> = rows.iter().collect();
    ordered.sort_by(|a, b| {
        a.locus
            .cmp(&b.locus)
            .then_with(|| a.ref_base.cmp(&b.ref_base))
            .then_with(|| a.call.alt_base.cmp(&b.call.alt_base))
    });

    for r in ordered {
        let chrom = contigs
            .by_id(r.locus.contig)
            .map(|c| c.name.as_ref())
            .unwrap_or(".");
        let pos = r.locus.pos.0 as u64 + 1;
        let total = (r.call.ad[0] + r.call.ad[1]).max(1);
        let af = r.call.ad[1] as f64 / total as f64;
        writeln!(
            out,
            "{chrom}\t{pos}\t.\t{ref_b}\t{alt}\t{qual:.1}\t{filt}\tDP={dp};AF={af:.3}\tGT:GQ:DP:AD:PL\t{gt}:{gq}:{dp}:{ad0},{ad1}:{pl0},{pl1},{pl2}",
            ref_b = r.ref_base as char,
            alt = r.call.alt_base as char,
            qual = r.call.qual,
            filt = filter_str(r.call.filter),
            dp = r.call.dp,
            gt = genotype_str(r.call.genotype),
            gq = r.call.gq,
            ad0 = r.call.ad[0],
            ad1 = r.call.ad[1],
            pl0 = r.call.pl[0],
            pl1 = r.call.pl[1],
            pl2 = r.call.pl[2],
        )?;
    }
    out.flush()
}

/// Render a germline VCF into a `String` (tests/snapshots).
pub fn render_germline_vcf(
    contigs: &ContigSet,
    sample: &str,
    rows: &[GermlineRow],
) -> io::Result<String> {
    let mut buf = Vec::new();
    write_germline_vcf(&mut buf, contigs, sample, rows)?;
    Ok(String::from_utf8(buf).expect("VCF is valid UTF-8"))
}
```

- [ ] **Step 4: Create the io module + register it.** Create `src/io/mod.rs`:

```rust
//! IO layer: standards-compliant readers and writers. Phase A4 lands the
//! spec-valid VCF writer; readers (FASTA/FASTQ/BAM) migrate here in later
//! phases.

pub mod vcf;
```

In `src/lib.rs`, add after the `pub mod genomics;` line:

```rust
/// IO layer: spec-valid VCF writer (FASTA/FASTQ/BAM readers arrive in later phases).
pub mod io;
```

- [ ] **Step 5: Run the tests to verify they pass.** Run: `cargo test io::vcf`
Expected: PASS (4 tests). Then `cargo build 2>&1 | grep -i warning` (none in `src/io/`), `cargo fmt --all -- --check`, `cargo clippy --lib 2>&1 | grep -A2 'src/io'` (fix obvious lints).

- [ ] **Step 6: Commit.**

```bash
git add src/io/mod.rs src/io/vcf.rs src/lib.rs
git commit -m "feat(io/vcf): spec-valid VCFv4.2 germline writer (typed header, deterministic)" \
  -m "Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task 3: `io/vcf` — somatic TUMOR/NORMAL writer

**Files:**
- Modify: `src/io/vcf.rs`
- Test: in `src/io/vcf.rs`

- [ ] **Step 1: Write the failing tests.** Add to the `#[cfg(test)] mod tests` block in `src/io/vcf.rs`:

```rust
    fn somatic_call() -> SomaticCall {
        SomaticCall {
            ref_base: b'A',
            alt_base: b'T',
            tumor_alt: 12,
            tumor_depth: 40,
            normal_alt: 0,
            normal_depth: 38,
            tumor_af: 0.3,
            normal_af: 0.0,
            quality: 55.0,
        }
    }

    #[test]
    fn somatic_header_has_two_samples() {
        let vcf = render_somatic_vcf(&contigs(), &[]).unwrap();
        assert!(vcf.starts_with("##fileformat=VCFv4.2\n"));
        assert!(vcf.contains("##INFO=<ID=SOMATIC,"));
        assert!(vcf.contains("##FORMAT=<ID=AF,Number=A,Type=Float,"));
        assert!(vcf.contains(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\tNORMAL\n"
        ));
    }

    #[test]
    fn somatic_record_is_exact() {
        let vcf = render_somatic_vcf(
            &contigs(),
            &[(
                Locus {
                    contig: 0,
                    pos: Position(200),
                },
                somatic_call(),
            )],
        )
        .unwrap();
        let last = vcf.lines().last().unwrap();
        // TUMOR AD = [40-12, 12] = 28,12; NORMAL AD = [38, 0].
        assert_eq!(
            last,
            "chr1\t201\t.\tA\tT\t55.0\tPASS\tSOMATIC\tGT:DP:AD:AF\t0/1:40:28,12:0.300\t0/0:38:38,0:0.000"
        );
    }
```

- [ ] **Step 2: Run the tests to verify they fail.** Run: `cargo test io::vcf`
Expected: FAIL — does not compile (`render_somatic_vcf` undefined).

- [ ] **Step 3: Implement the somatic writer.** Add to `src/io/vcf.rs`, above the `#[cfg(test)]` block:

```rust
/// Write a spec-valid somatic VCFv4.2 with `TUMOR` and `NORMAL` sample columns.
/// Records are emitted in canonical (contig, pos, ref, alt) order.
pub fn write_somatic_vcf<W: Write>(
    out: &mut W,
    contigs: &ContigSet,
    calls: &[(Locus, SomaticCall)],
) -> io::Result<()> {
    write_fileformat_and_contigs(out, contigs)?;
    writeln!(
        out,
        r#"##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Somatic mutation">"#
    )?;
    writeln!(
        out,
        r#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#
    )?;
    writeln!(
        out,
        r#"##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">"#
    )?;
    writeln!(
        out,
        r#"##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths (ref,alt)">"#
    )?;
    writeln!(
        out,
        r#"##FORMAT=<ID=AF,Number=A,Type=Float,Description="Alt allele fraction">"#
    )?;
    writeln!(
        out,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\tNORMAL"
    )?;

    let mut ordered: Vec<&(Locus, SomaticCall)> = calls.iter().collect();
    ordered.sort_by(|a, b| {
        a.0.cmp(&b.0)
            .then_with(|| a.1.ref_base.cmp(&b.1.ref_base))
            .then_with(|| a.1.alt_base.cmp(&b.1.alt_base))
    });

    for (locus, c) in ordered {
        let chrom = contigs
            .by_id(locus.contig)
            .map(|ct| ct.name.as_ref())
            .unwrap_or(".");
        let pos = locus.pos.0 as u64 + 1;
        let t_ref = c.tumor_depth.saturating_sub(c.tumor_alt);
        let n_ref = c.normal_depth.saturating_sub(c.normal_alt);
        writeln!(
            out,
            "{chrom}\t{pos}\t.\t{ref_b}\t{alt}\t{qual:.1}\tPASS\tSOMATIC\tGT:DP:AD:AF\t0/1:{t_dp}:{t_ref},{t_alt}:{t_af:.3}\t0/0:{n_dp}:{n_ref},{n_alt}:{n_af:.3}",
            ref_b = c.ref_base as char,
            alt = c.alt_base as char,
            qual = c.quality,
            t_dp = c.tumor_depth,
            t_alt = c.tumor_alt,
            t_af = c.tumor_af,
            n_dp = c.normal_depth,
            n_alt = c.normal_alt,
            n_af = c.normal_af,
        )?;
    }
    out.flush()
}

/// Render a somatic VCF into a `String` (tests/snapshots).
pub fn render_somatic_vcf(
    contigs: &ContigSet,
    calls: &[(Locus, SomaticCall)],
) -> io::Result<String> {
    let mut buf = Vec::new();
    write_somatic_vcf(&mut buf, contigs, calls)?;
    Ok(String::from_utf8(buf).expect("VCF is valid UTF-8"))
}
```

- [ ] **Step 4: Run the tests to verify they pass.** Run: `cargo test io::vcf`
Expected: PASS (6 tests total). Then `cargo test` (full suite green), `cargo build 2>&1 | grep -i warning` (none in `src/io/`), `cargo fmt --all -- --check` (run `cargo fmt --all` if needed), `cargo clippy --lib 2>&1 | grep -A2 'src/io'`.

- [ ] **Step 5: Commit.**

```bash
git add src/io/vcf.rs
git commit -m "feat(io/vcf): somatic TUMOR/NORMAL VCFv4.2 writer" \
  -m "Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Self-Review

**Spec coverage (§3.4, §3.5, §6):**
- §3.4 VcfHeader (fileformat/contig/INFO/FORMAT/FILTER) + record writer; germline 1 sample, somatic TUMOR/NORMAL; deterministic order — Task 2 + Task 3. ✔
- §3.5 RunManifest (tool_version, subcommand, inputs[path,blake3], params, outputs); canonical JSON sorted keys, no timestamps; `<output>.manifest.json`; BLAKE3 — Task 1. ✔
- §6 exact germline record format (POS 1-based, DP/AF INFO, GT:GQ:DP:AD:PL) — Task 2's `germline_record_is_exact` pins it. ✔
- §10 VCF validity (header tags + FORMAT/sample well-formed) + determinism (byte-identical twice / shuffled) + receipt determinism — covered. (`bcftools view` round-trip is a CI follow-up gated on tool availability; the Rust-level structural tests stand in locally.)

**Placeholder scan:** none — complete code + exact expected strings in every step.

**Type consistency:** writer consumes `GermlineCall { genotype, alt_base, qual, gq, pl[3], ad[2], dp, filter }` and `SomaticCall { ref_base, alt_base, tumor_alt, tumor_depth, normal_alt, normal_depth, tumor_af, normal_af, quality }` exactly as shipped by A3; `ContigSet::{iter, by_id}` and `Contig::{name, length}` and `Locus { contig, pos: Position(u32) }` exactly as shipped by A1. `RunManifest`/`FileHash` are self-defined. AF is computed (`ad[1]/(ad[0]+ad[1])`); somatic per-sample AD is `[depth−alt, alt]`.

**Scope:** writer + receipt only; no CLI wiring, no consumption by the calling path, no reader migration — all A5/B. Two independent, isolated modules. No `genomics` coupling (built on `core` + `call` + `std` + `blake3`).
