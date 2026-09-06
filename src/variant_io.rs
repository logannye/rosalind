//! Sequential, checked VCF/VCF.gz/BCF I/O for evidence selections and annotations.
//!
//! Inputs need no index and retain their record and allele order. HTSlib performs
//! native decoding before the declared header/record envelope is checked; these
//! are cooperative limits, not a guarantee against a transient native allocation.
//! Writers must be explicitly finished before an enclosing atomic artifact is
//! committed. Dropping a writer never reports successful completion.

use crate::core::ContigSet;
use crate::evidence::{EvidenceError, SnvSite};
use rust_htslib::bcf::{self, header::HeaderView, Read};
use rust_htslib::htslib;
use std::ffi::{CStr, CString};
use std::mem::size_of;
use std::path::Path;
use std::ptr::NonNull;

/// Declared per-header and per-decoded-record cooperative envelopes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct VariantLimits {
    /// Maximum formatted header bytes, including native formatting capacity.
    pub max_header_bytes: usize,
    /// Maximum decoded record capacities and, for VCF output, formatted capacity.
    pub max_record_bytes: usize,
}

impl Default for VariantLimits {
    fn default() -> Self {
        Self {
            max_header_bytes: 8 * 1024 * 1024,
            max_record_bytes: 1024 * 1024,
        }
    }
}

impl VariantLimits {
    fn validate(self) -> Result<(), EvidenceError> {
        if self.max_header_bytes == 0 || self.max_record_bytes == 0 {
            return Err(EvidenceError::InvalidRequest(
                "variant header and record envelopes must be positive".into(),
            ));
        }
        Ok(())
    }
}

/// A sequential reader with one reusable decoded record and checked parser flags.
#[derive(Debug)]
pub struct VariantReader {
    reader: bcf::Reader,
    record: bcf::Record,
    limits: VariantLimits,
    record_number: u64,
    failed: bool,
}

impl VariantReader {
    /// Open a local VCF, compressed VCF, or BCF without requiring an index.
    pub fn open(path: impl AsRef<Path>, limits: VariantLimits) -> Result<Self, EvidenceError> {
        limits.validate()?;
        check_input_file(path.as_ref())?;
        let reader = bcf::Reader::from_path(path.as_ref())
            .map_err(|e| invalid(format!("cannot open variant input: {e}")))?;
        // rust-htslib 0.44.1 does not check bcf_hdr_read's null result. Its
        // HeaderView destructor accepts null, but all other methods require it.
        if reader.header().inner.is_null() {
            return Err(invalid("cannot decode variant header"));
        }
        check_header(reader.header().inner, limits)?;
        let record = reader.empty_record();
        Ok(Self {
            reader,
            record,
            limits,
            record_number: 0,
            failed: false,
        })
    }

    /// The unchanged input header, including sample and INFO/FORMAT definitions.
    pub fn header(&self) -> &HeaderView {
        self.reader.header()
    }

    /// Number of records decoded so far, including a record that failed validation.
    pub fn record_number(&self) -> u64 {
        self.record_number
    }

    /// Read the next record in original order, reusing the same native allocation.
    /// Any parse or envelope failure poisons the reader rather than allowing a
    /// caller to continue with header definitions invented by HTSlib recovery.
    pub fn read(&mut self) -> Result<Option<&mut bcf::Record>, EvidenceError> {
        if self.failed {
            return Err(invalid(
                "variant reader cannot continue after a failed record",
            ));
        }
        let result = match self.reader.read(&mut self.record) {
            None => return Ok(None),
            Some(result) => result,
        };
        self.record_number = self
            .record_number
            .checked_add(1)
            .ok_or(EvidenceError::CounterOverflow)?;
        let checked = result
            .map_err(|e| invalid(format!("variant record {}: {e}", self.record_number)))
            .and_then(|()| {
                // SAFETY: rust-htslib initialized this record. Its read method
                // ignores the unpack return value, so check it explicitly too.
                if unsafe { htslib::bcf_unpack(self.record.inner, htslib::BCF_UN_ALL as i32) } < 0 {
                    return Err(invalid("cannot unpack variant input record"));
                }
                check_record(self.record.inner(), self.limits)?;
                check_record_header(&self.record)
            });
        if let Err(error) = checked {
            self.failed = true;
            return Err(error);
        }
        Ok(Some(&mut self.record))
    }
}

/// Resolve a single typed record to an A/C/G/T SNV, preserving its ALT order.
/// Selection normalization subsequently unions repeated loci and ALT alleles.
pub fn parse_snv_record(
    record: &bcf::Record,
    contigs: &ContigSet,
) -> Result<SnvSite, EvidenceError> {
    if record.inner().errcode != 0 {
        return Err(invalid(
            "variant record has unresolved HTSlib parser errors",
        ));
    }
    let rid = record
        .rid()
        .ok_or_else(|| invalid("variant has no contig"))?;
    let name = record
        .header()
        .rid2name(rid)
        .map_err(|e| invalid(format!("invalid variant contig: {e}")))?;
    let name = std::str::from_utf8(name).map_err(|_| invalid("non-UTF-8 variant contig"))?;
    let contig = contigs
        .by_name(name)
        .ok_or_else(|| invalid(format!("unknown variant contig {name}")))?;
    let position = u32::try_from(record.pos()).map_err(|_| invalid("POS outside reference"))?;
    if position >= contig.length {
        return Err(invalid("POS outside reference"));
    }
    let alleles = record.alleles();
    if alleles.len() < 2 {
        return Err(invalid("SNV record requires at least one ALT"));
    }
    let base = |allele: &[u8]| {
        if allele.len() != 1 || !b"ACGT".contains(&allele[0].to_ascii_uppercase()) {
            return Err(invalid(
                "only A/C/G/T SNV REF and ALT alleles are supported",
            ));
        }
        Ok(allele[0].to_ascii_uppercase())
    };
    let reference = base(alleles[0])?;
    if alleles.len() > 4 {
        return Err(invalid(
            "SNV record cannot have more than three distinct ALT bases",
        ));
    }
    let mut seen = [false; 256];
    let alternates = alleles[1..]
        .iter()
        .map(|allele| {
            let alternate = base(allele)?;
            if alternate == reference {
                return Err(invalid("ALT equals REF"));
            }
            if std::mem::replace(&mut seen[alternate as usize], true) {
                return Err(invalid("duplicate ALT allele within a variant record"));
            }
            Ok(alternate)
        })
        .collect::<Result<Vec<_>, _>>()?;
    Ok(SnvSite {
        contig: contig.id,
        position,
        reference,
        alternates,
    })
}

/// Explicit output encoding, independent of a temporary artifact's suffix.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VariantFormat {
    /// Uncompressed VCF text.
    Vcf,
    /// BGZF-compressed VCF text.
    VcfGz,
    /// BGZF-compressed BCF.
    Bcf,
}

/// Checked native writer with a reusable copy for header translation/annotations.
///
/// The original Rust record remains associated with its input header. Account for
/// input/output/translation headers, native I/O, the input record, the copied
/// record, and VCF formatting scratch when planning an annotation process.
#[derive(Debug)]
pub struct CheckedVariantWriter {
    file: Option<NonNull<htslib::htsFile>>,
    header: HeaderView,
    source: Option<SourceHeader>,
    scratch: NativeRecord,
    text: NativeString,
    format: VariantFormat,
    limits: VariantLimits,
    failed: bool,
}

impl CheckedVariantWriter {
    /// Create an output and check header serialization/write errors immediately.
    pub fn create(
        path: impl AsRef<Path>,
        header: &bcf::Header,
        format: VariantFormat,
        limits: VariantLimits,
    ) -> Result<Self, EvidenceError> {
        limits.validate()?;
        let path = checked_path(path.as_ref())?;
        if header.inner.is_null() {
            return Err(invalid("null variant output header"));
        }
        check_header(header.inner, limits)?;
        // SAFETY: the borrowed Header owns a live bcf_hdr_t; the duplicate becomes
        // exclusively owned by HeaderView and is destroyed exactly once.
        let cloned = unsafe { htslib::bcf_hdr_dup(header.inner) };
        if cloned.is_null() {
            return Err(io_error("cannot duplicate variant output header"));
        }
        let header = HeaderView::new(cloned);
        check_header(header.inner, limits)?;
        let mode = match format {
            VariantFormat::Vcf => c"w",
            VariantFormat::VcfGz => c"wz",
            VariantFormat::Bcf => c"wb",
        };
        let scratch = NativeRecord::new()?;
        // SAFETY: both C strings are NUL terminated and live for the call.
        let file = NonNull::new(unsafe { htslib::hts_open(path.as_ptr(), mode.as_ptr()) })
            .ok_or_else(|| io_error("cannot open variant output"))?;
        let mut writer = Self {
            file: Some(file),
            header,
            source: None,
            scratch,
            text: NativeString::default(),
            format,
            limits,
            failed: false,
        };
        // SAFETY: file and header are live, uniquely owned native objects.
        if unsafe { htslib::bcf_hdr_write(file.as_ptr(), writer.header.inner) } < 0 {
            writer.failed = true;
            return Err(io_error("cannot write variant output header"));
        }
        Ok(writer)
    }

    /// The owned output header, including added annotation definitions.
    pub fn header(&self) -> &HeaderView {
        &self.header
    }

    /// Write a translated copy without changing the original record or its header.
    pub fn write(&mut self, record: &bcf::Record) -> Result<(), EvidenceError> {
        self.write_with_info(record, &[], &[])
    }

    /// Add typed INFO fields to a translated copy and check the complete write.
    /// Keys must already exist with the matching type in the output header.
    /// String values are comma-separated VCF elements and must not contain
    /// delimiters, whitespace, or NUL. Integer narrowing is the caller's duty.
    pub fn write_with_info(
        &mut self,
        record: &bcf::Record,
        integers: &[(&[u8], &[i32])],
        strings: &[(&[u8], &[&[u8]])],
    ) -> Result<(), EvidenceError> {
        if self.failed {
            return Err(io_error(
                "variant writer cannot continue after a failed write",
            ));
        }
        let result = self.write_record(record, integers, strings);
        if result.is_err() {
            self.failed = true;
        }
        result
    }

    fn write_record(
        &mut self,
        record: &bcf::Record,
        integers: &[(&[u8], &[i32])],
        strings: &[(&[u8], &[&[u8]])],
    ) -> Result<(), EvidenceError> {
        let source_bytes = check_record(record.inner(), self.limits)?;
        check_record_header(record)?;
        self.bind_source(record)?;
        let mut annotation_bytes = 0usize;
        for (key, values) in integers {
            self.check_info_type(key, bcf::header::TagType::Integer)?;
            annotation_bytes = annotation_bytes
                .checked_add(
                    values
                        .len()
                        .checked_mul(size_of::<i32>())
                        .ok_or(EvidenceError::CounterOverflow)?,
                )
                .ok_or(EvidenceError::CounterOverflow)?;
        }
        for (key, values) in strings {
            self.check_info_type(key, bcf::header::TagType::String)?;
            for value in *values {
                if value.is_empty()
                    || value
                        .iter()
                        .any(|b| b.is_ascii_whitespace() || b"\0,;=".contains(b))
                {
                    return Err(invalid("invalid VCF INFO string element"));
                }
                annotation_bytes = annotation_bytes
                    .checked_add(
                        value
                            .len()
                            .checked_add(1)
                            .ok_or(EvidenceError::CounterOverflow)?,
                    )
                    .ok_or(EvidenceError::CounterOverflow)?;
            }
        }
        check_limit(
            source_bytes
                .checked_add(annotation_bytes)
                .ok_or(EvidenceError::CounterOverflow)?,
            self.limits.max_record_bytes,
            "record plus annotations",
        )?;
        // SAFETY: source, reusable scratch, and both headers are live. bcf_copy
        // synchronizes native packed storage but preserves source scientific data;
        // translation/INFO updates apply only to the exclusive scratch record.
        unsafe {
            if htslib::bcf_copy(self.scratch.0.as_ptr(), record.inner).is_null() {
                return Err(io_error("cannot copy variant output record"));
            }
            if htslib::bcf_translate(
                self.header.inner,
                self.source.as_ref().expect("bound source").translated.inner,
                self.scratch.0.as_ptr(),
            ) < 0
            {
                return Err(invalid("cannot translate variant record to output header"));
            }
        }
        for (key, values) in integers {
            let key = checked_key(key)?;
            self.update_info(
                &key,
                values.as_ptr().cast(),
                values.len(),
                htslib::BCF_HT_INT,
            )?;
        }
        for (key, values) in strings {
            let key = checked_key(key)?;
            // BCF_HT_STR uses strlen: its count is not a byte-length bound.
            // Include the NUL already reserved by annotation_bytes above.
            let joined = CString::new(values.join(&b','))
                .map_err(|_| invalid("variant INFO string contains NUL"))?;
            self.update_info(
                &key,
                joined.as_ptr().cast(),
                usize::from(!values.is_empty()),
                htslib::BCF_HT_STR,
            )?;
        }
        // SAFETY: the scratch record remains owned and initialized. Unpacking
        // makes every capacity visible to the envelope accounting below.
        if unsafe { htslib::bcf_unpack(self.scratch.0.as_ptr(), htslib::BCF_UN_ALL as i32) } < 0 {
            return Err(invalid("cannot unpack annotated variant record"));
        }
        check_record(self.scratch.get(), self.limits)?;
        if self.format != VariantFormat::Bcf {
            self.text.0.l = 0;
            // SAFETY: header/record/string are live; native kstring owns its buffer.
            if unsafe {
                htslib::vcf_format(self.header.inner, self.scratch.0.as_ptr(), &mut self.text.0)
            } < 0
            {
                return Err(invalid("cannot format annotated VCF record"));
            }
            check_limit(
                usize::try_from(self.text.0.m).map_err(|_| EvidenceError::CounterOverflow)?,
                self.limits.max_record_bytes,
                "formatted VCF record",
            )?;
        }
        // SAFETY: the file is still open and all translated IDs use its header.
        if unsafe {
            htslib::bcf_write(
                self.file.expect("open writer").as_ptr(),
                self.header.inner,
                self.scratch.0.as_ptr(),
            )
        } < 0
        {
            return Err(io_error("cannot write variant output record"));
        }
        // BCF serialization can grow packed buffers after INFO updates.
        check_record(self.scratch.get(), self.limits)?;
        Ok(())
    }

    fn bind_source(&mut self, record: &bcf::Record) -> Result<(), EvidenceError> {
        if self
            .source
            .as_ref()
            .is_some_and(|source| source.anchor.header().inner == record.header().inner)
        {
            return Ok(());
        }
        check_header(record.header().inner, self.limits)?;
        check_output_header(record.header(), &self.header)?;
        // SAFETY: duplicate a live header. Translation tables belong to this
        // writer; bcf_translate must never cache them in the shared input header.
        let cloned = unsafe { htslib::bcf_hdr_dup(record.header().inner) };
        if cloned.is_null() {
            return Err(io_error("cannot duplicate variant translation header"));
        }
        let translated = HeaderView::new(cloned);
        let empty = NativeRecord::new()?;
        let mut anchor = record.clone();
        // Keep the source Rc<HeaderView> alive without retaining a record payload.
        // SAFETY: the new empty record replaces one owned native record; the old
        // allocation is freed once, and anchor's destructor owns the replacement.
        unsafe {
            let old = std::mem::replace(&mut anchor.inner, empty.0.as_ptr());
            std::mem::forget(empty);
            htslib::bcf_destroy(old);
        }
        self.source = Some(SourceHeader { anchor, translated });
        Ok(())
    }

    fn check_info_type(
        &self,
        key: &[u8],
        expected: bcf::header::TagType,
    ) -> Result<(), EvidenceError> {
        checked_key(key)?;
        let (actual, _) = self
            .header
            .info_type(key)
            .map_err(|e| invalid(format!("undefined output INFO key: {e}")))?;
        if actual != expected {
            return Err(invalid("output INFO key has incompatible type"));
        }
        Ok(())
    }

    fn update_info(
        &mut self,
        key: &CString,
        values: *const libc::c_void,
        count: usize,
        kind: u32,
    ) -> Result<(), EvidenceError> {
        let count = i32::try_from(count).map_err(|_| EvidenceError::CounterOverflow)?;
        // SAFETY: the caller's typed slice or NUL-terminated CString lives for
        // the call. HTSlib copies it into the exclusive scratch record.
        if unsafe {
            htslib::bcf_update_info(
                self.header.inner,
                self.scratch.0.as_ptr(),
                key.as_ptr(),
                values,
                count,
                kind as i32,
            )
        } < 0
        {
            return Err(invalid("cannot update variant INFO value"));
        }
        Ok(())
    }

    /// Flush and close explicitly. A prior write failure always prevents success.
    pub fn finish(mut self) -> Result<(), EvidenceError> {
        let file = self.file.take().expect("open writer");
        // SAFETY: take transfers the sole handle, preventing Drop from reclosing.
        let result = unsafe { htslib::hts_close(file.as_ptr()) };
        if result < 0 {
            return Err(io_error("cannot flush/close variant output"));
        }
        if self.failed {
            return Err(io_error("variant output had an earlier failed write"));
        }
        Ok(())
    }
}

impl Drop for CheckedVariantWriter {
    fn drop(&mut self) {
        if let Some(file) = self.file.take() {
            // SAFETY: best-effort abandonment closes the sole remaining handle.
            // Only explicit finish can report successful artifact completion.
            unsafe { htslib::hts_close(file.as_ptr()) };
        }
    }
}

#[derive(Debug)]
struct SourceHeader {
    // An empty Record retains the original private Rc<HeaderView>, preventing a
    // later reader from reusing its address while this translation cache is live.
    anchor: bcf::Record,
    translated: HeaderView,
}

#[derive(Debug)]
struct NativeRecord(NonNull<htslib::bcf1_t>);
impl NativeRecord {
    fn new() -> Result<Self, EvidenceError> {
        // SAFETY: allocation has no arguments; ownership is transferred on success.
        NonNull::new(unsafe { htslib::bcf_init() })
            .map(Self)
            .ok_or_else(|| io_error("cannot allocate variant record"))
    }
    fn get(&self) -> &htslib::bcf1_t {
        // SAFETY: this wrapper exclusively owns an initialized live record.
        unsafe { self.0.as_ref() }
    }
}
impl Drop for NativeRecord {
    fn drop(&mut self) {
        // SAFETY: this is the sole owning wrapper.
        unsafe { htslib::bcf_destroy(self.0.as_ptr()) };
    }
}

#[derive(Debug)]
struct NativeString(htslib::kstring_t);
impl Default for NativeString {
    fn default() -> Self {
        Self(htslib::kstring_t {
            l: 0,
            m: 0,
            s: std::ptr::null_mut(),
        })
    }
}
impl Drop for NativeString {
    fn drop(&mut self) {
        // SAFETY: HTSlib kstrings use malloc/realloc; free also accepts null.
        unsafe { libc::free(self.0.s.cast()) };
    }
}

fn check_header(
    header: *mut htslib::bcf_hdr_t,
    limits: VariantLimits,
) -> Result<(), EvidenceError> {
    let mut text = NativeString::default();
    // SAFETY: synchronizing cached IDs before duplication also captures samples
    // just appended to a Header. HTSlib's formatter alone does not synchronize.
    if unsafe { htslib::bcf_hdr_sync(header) } < 0 {
        return Err(invalid("cannot synchronize variant header"));
    }
    // SAFETY: callers validate nonnull, live header; text owns a native kstring.
    if unsafe { htslib::bcf_hdr_format(header, 0, &mut text.0) } < 0 {
        return Err(invalid("cannot format variant header"));
    }
    check_limit(
        usize::try_from(text.0.m).map_err(|_| EvidenceError::CounterOverflow)?,
        limits.max_header_bytes,
        "variant header",
    )?;
    if text.0.l != 0 {
        let length = usize::try_from(text.0.l).map_err(|_| EvidenceError::CounterOverflow)?;
        // SAFETY: successful formatting initialized exactly l bytes in kstring.
        let bytes = unsafe { std::slice::from_raw_parts(text.0.s.cast::<u8>(), length) };
        std::str::from_utf8(bytes).map_err(|_| invalid("variant header is not UTF-8"))?;
    }
    Ok(())
}

fn check_output_header(source: &HeaderView, output: &HeaderView) -> Result<(), EvidenceError> {
    let sample_count = source.sample_count();
    // rust-htslib 0.44.1's samples() constructs a slice from a null pointer for
    // some valid zero-sample headers; avoid calling it for that case.
    if sample_count != output.sample_count()
        || (sample_count != 0 && source.samples() != output.samples())
    {
        return Err(invalid(
            "variant output must preserve sample names and order",
        ));
    }
    for rid in 0..source.contig_count() {
        let name = source.rid2name(rid).map_err(|e| invalid(e.to_string()))?;
        let output_rid = output
            .name2rid(name)
            .map_err(|_| invalid("variant output is missing an input contig"))?;
        // SAFETY: name/rid lookups above validate the live contig dictionary
        // entries. HTSlib stores declared contig length in idinfo.info[0].
        let same_length = unsafe {
            let source_value =
                (*(*source.inner).id[htslib::BCF_DT_CTG as usize].add(rid as usize)).val;
            let output_value =
                (*(*output.inner).id[htslib::BCF_DT_CTG as usize].add(output_rid as usize)).val;
            !source_value.is_null()
                && !output_value.is_null()
                && (*source_value).info[0] == (*output_value).info[0]
        };
        if !same_length {
            return Err(invalid("variant output changes an input contig length"));
        }
    }
    // SAFETY: both headers are initialized and validated as UTF-8. Enumerate
    // HTSlib's live dictionary entries; removed entries have a null key.
    unsafe {
        let source_raw = &*source.inner;
        for index in 0..source_raw.n[htslib::BCF_DT_ID as usize] as usize {
            let pair = &*source_raw.id[htslib::BCF_DT_ID as usize].add(index);
            if pair.key.is_null() {
                continue;
            }
            let key = CStr::from_ptr(pair.key).to_bytes();
            output
                .name_to_id(key)
                .map_err(|_| invalid("variant output is missing an input tag"))?;
            for (source_type, output_type) in [
                (source.info_type(key), output.info_type(key)),
                (source.format_type(key), output.format_type(key)),
            ] {
                if let Ok(expected) = source_type {
                    if output_type.ok() != Some(expected) {
                        return Err(invalid(
                            "variant output changes an input INFO/FORMAT definition",
                        ));
                    }
                }
            }
            let source_filter = htslib::bcf_hdr_get_hrec(
                source.inner,
                htslib::BCF_HL_FLT as i32,
                c"ID".as_ptr(),
                pair.key,
                std::ptr::null(),
            );
            if !source_filter.is_null()
                && htslib::bcf_hdr_get_hrec(
                    output.inner,
                    htslib::BCF_HL_FLT as i32,
                    c"ID".as_ptr(),
                    pair.key,
                    std::ptr::null(),
                )
                .is_null()
            {
                return Err(invalid(
                    "variant output is missing an input FILTER definition",
                ));
            }
        }
    }
    Ok(())
}

fn check_record_header(record: &bcf::Record) -> Result<(), EvidenceError> {
    if record.inner().n_sample() != record.header().sample_count() {
        return Err(invalid(
            "variant record sample count does not match its header",
        ));
    }
    if record
        .rid()
        .is_none_or(|rid| rid >= record.header().contig_count())
    {
        return Err(invalid("variant record contig ID is outside its header"));
    }
    Ok(())
}

fn check_record(record: &htslib::bcf1_t, limits: VariantLimits) -> Result<usize, EvidenceError> {
    // Bindgen names this pointee differently on Linux and macOS. Infer its
    // native layout from the field type without dereferencing the pointer.
    fn pointee_size<T>(_: *const T) -> usize {
        size_of::<T>()
    }
    if record.errcode != 0 {
        return Err(invalid(format!(
            "variant parser error flags {} (undefined contig/tag or malformed record)",
            record.errcode
        )));
    }
    let d = &record.d;
    let mut bytes = size_of::<htslib::bcf1_t>();
    for capacity in [record.shared.m, record.indiv.m] {
        bytes = bytes
            .checked_add(usize::try_from(capacity).map_err(|_| EvidenceError::CounterOverflow)?)
            .ok_or(EvidenceError::CounterOverflow)?;
    }
    for (count, width) in [
        (d.m_id, 1),
        (d.m_als, 1),
        (d.m_allele, size_of::<*mut libc::c_char>()),
        (d.m_flt, size_of::<libc::c_int>()),
        (d.m_info, size_of::<htslib::bcf_info_t>()),
        (d.m_fmt, size_of::<htslib::bcf_fmt_t>()),
        (d.n_var, pointee_size(d.var)),
    ] {
        let count =
            usize::try_from(count).map_err(|_| invalid("negative variant allocation capacity"))?;
        bytes = bytes
            .checked_add(
                count
                    .checked_mul(width)
                    .ok_or(EvidenceError::CounterOverflow)?,
            )
            .ok_or(EvidenceError::CounterOverflow)?;
    }
    // Updated INFO/FORMAT values can own independent buffers outside shared/indiv.
    check_limit(bytes, limits.max_record_bytes, "decoded variant record")?;
    let info_count = if record.unpacked & htslib::BCF_UN_INFO as i32 != 0 {
        record.n_info() as usize
    } else {
        0
    };
    let format_count = if record.unpacked & htslib::BCF_UN_FMT as i32 != 0 {
        record.n_fmt() as usize
    } else {
        0
    };
    if (info_count > 0 && (d.info.is_null() || info_count > d.m_info as usize))
        || (format_count > 0 && (d.fmt.is_null() || format_count > d.m_fmt as usize))
    {
        return Err(invalid("variant decoded field array is inconsistent"));
    }
    // SAFETY: capacities/pointers come from an unpacked live native record. Only
    // initialized entries (n_info/n_fmt), rather than spare capacity, are read.
    unsafe {
        for i in 0..info_count {
            let info = &*d.info.add(i);
            if info.vptr_free() != 0 {
                bytes = bytes
                    .checked_add(info.vptr_len as usize + info.vptr_off() as usize)
                    .ok_or(EvidenceError::CounterOverflow)?;
            }
        }
        for i in 0..format_count {
            let format = &*d.fmt.add(i);
            if format.p_free() != 0 {
                bytes = bytes
                    .checked_add(format.p_len as usize + format.p_off() as usize)
                    .ok_or(EvidenceError::CounterOverflow)?;
            }
        }
    }
    check_limit(bytes, limits.max_record_bytes, "decoded variant record")?;
    Ok(bytes)
}

fn checked_path(path: &Path) -> Result<CString, EvidenceError> {
    let text = path
        .to_str()
        .ok_or_else(|| invalid("variant path is not UTF-8"))?;
    if text == "-" {
        return Err(invalid(
            "variant artifacts require a local path, not standard I/O",
        ));
    }
    CString::new(text).map_err(|_| invalid("variant path contains NUL"))
}

fn check_input_file(path: &Path) -> Result<(), EvidenceError> {
    let cpath = checked_path(path)?;
    if !std::fs::metadata(path)?.is_file() {
        return Err(invalid("variant input must be a regular local file"));
    }
    // The pinned high-level Reader hides its htsFile handle. A short independent
    // open checks the BGZF end marker before sequential decoding, so a truncated
    // stream cannot be mistaken for normal EOF by the high-level iterator.
    // SAFETY: path/mode are live C strings; this scope owns the returned handle.
    let file = NonNull::new(unsafe { htslib::hts_open(cpath.as_ptr(), c"r".as_ptr()) })
        .ok_or_else(|| invalid("cannot inspect variant input"))?;
    // SAFETY: both calls use the live handle, closed exactly once here.
    let eof = unsafe { htslib::hts_check_EOF(file.as_ptr()) };
    let closed = unsafe { htslib::hts_close(file.as_ptr()) };
    if !matches!(eof, 1 | 3) {
        return Err(invalid(
            "variant input has no valid BGZF end marker or cannot be checked",
        ));
    }
    if closed < 0 {
        return Err(io_error("cannot close variant input preflight handle"));
    }
    Ok(())
}

fn checked_key(key: &[u8]) -> Result<CString, EvidenceError> {
    if key.is_empty() {
        return Err(invalid("empty variant INFO key"));
    }
    CString::new(key).map_err(|_| invalid("variant INFO key contains NUL"))
}

fn check_limit(actual: usize, maximum: usize, label: &str) -> Result<(), EvidenceError> {
    if actual > maximum {
        return Err(EvidenceError::RecordLimit(format!(
            "{label} requires {actual} bytes; declared {maximum}"
        )));
    }
    Ok(())
}

fn invalid(message: impl Into<String>) -> EvidenceError {
    EvidenceError::InvalidInput(message.into())
}

fn io_error(message: &str) -> EvidenceError {
    EvidenceError::Io(std::io::Error::other(message))
}
