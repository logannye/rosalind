//! Minimal, read-only memory mapping helper.
//!
//! We avoid adding new dependencies for mmap (e.g. `memmap2`) and instead use
//! `libc::mmap` directly. This is intentionally small and only supports the
//! needs of read-only reference/index files.

use std::fs::File;
use std::io;
use std::os::fd::AsRawFd;
use std::ptr::NonNull;

/// A read-only memory-mapped file.
///
/// # Safety invariants
/// - `ptr..ptr+len` is a valid read-only mapping for the lifetime of `Self`.
/// - The mapping is unmapped on drop.
#[derive(Debug)]
pub struct MmapReadOnly {
    ptr: NonNull<u8>,
    len: usize,
}

impl MmapReadOnly {
    /// Memory-map the entire file as read-only.
    pub fn map(file: &File) -> io::Result<Self> {
        let len = file.metadata()?.len() as usize;
        if len == 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "cannot mmap empty file",
            ));
        }

        // SAFETY: `mmap` is called with a valid FD, length > 0, and read-only flags.
        // We check for MAP_FAILED and wrap the pointer in NonNull on success.
        unsafe {
            let addr = libc::mmap(
                std::ptr::null_mut(),
                len,
                libc::PROT_READ,
                libc::MAP_PRIVATE,
                file.as_raw_fd(),
                0,
            );
            if addr == libc::MAP_FAILED {
                return Err(io::Error::last_os_error());
            }
            let ptr = NonNull::new(addr as *mut u8).ok_or_else(|| {
                io::Error::new(io::ErrorKind::Other, "mmap returned null pointer")
            })?;
            Ok(Self { ptr, len })
        }
    }

    /// Expose the mapped bytes.
    pub fn as_bytes(&self) -> &[u8] {
        // SAFETY: the mapping is valid for `len` bytes for the lifetime of `Self`.
        unsafe { std::slice::from_raw_parts(self.ptr.as_ptr(), self.len) }
    }
}

impl Drop for MmapReadOnly {
    fn drop(&mut self) {
        // SAFETY: `ptr` was created by `mmap` and is valid for `len` bytes.
        unsafe {
            libc::munmap(self.ptr.as_ptr() as *mut libc::c_void, self.len);
        }
    }
}


