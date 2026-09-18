//! Narrow indexed alignment reader with explicit native ownership ordering.
//!
//! rust-htslib 0.44.1's IndexedReader closes its CRAM handle before its index
//! is destroyed. CRAM indexes borrow that handle, so their destructor then uses
//! freed memory (upstream rust-bio/rust-htslib#518). Keep the pinned decoder ABI
//! and use the existing sequential Reader to own the file, but destroy our
//! iterator and index before dropping that Reader, including on error paths.

use rust_htslib::{bam, bam::Read, errors::Error, htslib};
use std::{ffi::CString, path::Path, ptr::NonNull};

#[derive(Debug)]
pub(super) struct IndexedAlignmentReader {
    reader: bam::Reader,
    index: NonNull<htslib::hts_idx_t>,
    iterator: Option<NonNull<htslib::hts_itr_t>>,
}

// SAFETY: all native pointers are exclusively owned, never exposed, and used
// only with mutable access. Moving the owner also moves its bam::Reader (Send).
// No native decoder thread is started. This preserves EvidenceEngine's existing
// Send contract while deliberately not making the reader Sync.
unsafe impl Send for IndexedAlignmentReader {}

impl IndexedAlignmentReader {
    pub(super) fn open(path: &Path, index_path: Option<&Path>) -> Result<Self, Error> {
        let path_string = native_path(path)?;
        let index_string = index_path.map(native_path).transpose()?;
        let reader = bam::Reader::from_path(path)?;
        // SAFETY: Reader owns a live htsFile. Both C strings remain live for the
        // call; HTSlib returns a separately owned index or null on failure.
        let index = unsafe {
            htslib::sam_index_load2(
                reader.htsfile(),
                path_string.as_ptr(),
                index_string
                    .as_ref()
                    .map_or(std::ptr::null(), |p| p.as_ptr()),
            )
        };
        let index = NonNull::new(index).ok_or_else(|| Error::BamInvalidIndex {
            target: path.display().to_string(),
        })?;
        Ok(Self {
            reader,
            index,
            iterator: None,
        })
    }

    pub(super) fn header(&self) -> &bam::HeaderView {
        self.reader.header()
    }

    pub(super) fn set_reference(&mut self, path: &Path) -> Result<(), Error> {
        self.reader.set_reference(path)
    }

    fn clear_iterator(&mut self) {
        if let Some(iterator) = self.iterator.take() {
            // SAFETY: this reader exclusively owns the iterator; take prevents
            // a second destruction, and its index and file are still live.
            unsafe { htslib::hts_itr_destroy(iterator.as_ptr()) };
        }
    }

    pub(super) fn fetch(&mut self, (tid, start, end): (u32, i64, i64)) -> Result<(), Error> {
        self.clear_iterator();
        let tid = i32::try_from(tid).map_err(|_| Error::Fetch)?;
        if start < 0 || end < start {
            return Err(Error::Fetch);
        }
        // SAFETY: the index and the file it may borrow remain owned by self.
        self.iterator =
            NonNull::new(unsafe { htslib::sam_itr_queryi(self.index.as_ptr(), tid, start, end) });
        self.iterator.map(|_| ()).ok_or(Error::Fetch)
    }

    pub(super) fn read(&mut self, record: &mut bam::Record) -> Option<Result<(), Error>> {
        let iterator = self.iterator?;
        let file = self.reader.htsfile();
        // SAFETY: self owns a live file, index and iterator; record owns its
        // writable BAM payload. This is the pinned rust-htslib iterator call,
        // including the htsFile context needed by CRAM. The evidence consumer
        // uses numeric tids and never requires a borrowed record header or a
        // cached CIGAR, so no additional header allocation is needed.
        let result = unsafe {
            htslib::hts_itr_next(
                (*file).fp.bgzf,
                iterator.as_ptr(),
                (record.inner_mut() as *mut htslib::bam1_t).cast(),
                file.cast(),
            )
        };
        match result {
            -1 => None,
            -2 => Some(Err(Error::BamTruncatedRecord)),
            result if result < 0 => Some(Err(Error::BamInvalidRecord)),
            _ => Some(Ok(())),
        }
    }
}

impl Drop for IndexedAlignmentReader {
    fn drop(&mut self) {
        self.clear_iterator();
        // SAFETY: the exclusively owned index must be destroyed while Reader's
        // CRAM file is live. Rust drops the reader field only after this method.
        unsafe { htslib::hts_idx_destroy(self.index.as_ptr()) };
    }
}

fn native_path(path: &Path) -> Result<CString, Error> {
    path.to_str()
        .and_then(|path| CString::new(path).ok())
        .ok_or_else(|| Error::BamInvalidIndex {
            target: path.display().to_string(),
        })
}
