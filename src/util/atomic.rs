//! Transactional file creation for user-facing artifacts and receipts.

use std::fs::{File, OpenOptions};
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

static TEMP_COUNTER: AtomicU64 = AtomicU64::new(0);

/// A same-directory temporary file that can be atomically committed.
#[derive(Debug)]
pub struct AtomicFile {
    destination: PathBuf,
    temporary: PathBuf,
    file: Option<File>,
    committed: bool,
}

impl AtomicFile {
    /// Create a unique sibling temporary file for `destination`.
    pub fn create(destination: &Path) -> io::Result<Self> {
        let parent = destination.parent().unwrap_or_else(|| Path::new("."));
        let name = destination
            .file_name()
            .and_then(|name| name.to_str())
            .unwrap_or("artifact");
        for _ in 0..128 {
            let counter = TEMP_COUNTER.fetch_add(1, Ordering::Relaxed);
            let temporary = parent.join(format!(
                ".{name}.rosalind-{}-{counter}.partial",
                std::process::id()
            ));
            match OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(&temporary)
            {
                Ok(file) => {
                    return Ok(Self {
                        destination: destination.to_path_buf(),
                        temporary,
                        file: Some(file),
                        committed: false,
                    })
                }
                Err(error) if error.kind() == io::ErrorKind::AlreadyExists => continue,
                Err(error) => return Err(error),
            }
        }
        Err(io::Error::new(
            io::ErrorKind::AlreadyExists,
            "could not reserve a unique transactional output path",
        ))
    }

    /// Borrow the writable temporary file.
    pub fn file_mut(&mut self) -> &mut File {
        self.file.as_mut().expect("atomic file is still open")
    }

    /// Path holding bytes before commit.
    pub fn temporary_path(&self) -> &Path {
        &self.temporary
    }

    /// Flush, sync, close, and atomically move the file to its destination.
    pub fn commit(mut self, replace: bool) -> io::Result<PathBuf> {
        let destination = self.destination.clone();
        self.finish(&destination, replace)
    }

    /// Commit the bytes to a different sibling path, used for governed partial output.
    pub fn commit_as(mut self, destination: &Path, replace: bool) -> io::Result<PathBuf> {
        self.finish(destination, replace)
    }

    fn finish(&mut self, destination: &Path, replace: bool) -> io::Result<PathBuf> {
        if let Some(mut file) = self.file.take() {
            file.flush()?;
            file.sync_all()?;
        }
        if !replace && destination.exists() {
            return Err(io::Error::new(
                io::ErrorKind::AlreadyExists,
                format!("destination already exists: {}", destination.display()),
            ));
        }
        #[cfg(windows)]
        if replace && destination.exists() {
            std::fs::remove_file(destination)?;
        }
        std::fs::rename(&self.temporary, destination)?;
        self.committed = true;
        Ok(destination.to_path_buf())
    }
}

impl Drop for AtomicFile {
    fn drop(&mut self) {
        if !self.committed {
            self.file.take();
            let _ = std::fs::remove_file(&self.temporary);
        }
    }
}

/// Atomically write a complete small file, refusing replacement unless requested.
pub fn write_atomic(path: &Path, bytes: &[u8], replace: bool) -> io::Result<()> {
    let mut file = AtomicFile::create(path)?;
    file.file_mut().write_all(bytes)?;
    file.commit(replace)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn path(name: &str) -> PathBuf {
        std::env::temp_dir().join(format!(
            "rosalind-atomic-{name}-{}-{}",
            std::process::id(),
            TEMP_COUNTER.fetch_add(1, Ordering::Relaxed)
        ))
    }

    #[test]
    fn commit_is_atomic_and_create_new_refuses_replacement() {
        let output = path("commit");
        write_atomic(&output, b"first", false).unwrap();
        let error = write_atomic(&output, b"second", false).unwrap_err();
        assert_eq!(error.kind(), io::ErrorKind::AlreadyExists);
        assert_eq!(std::fs::read(&output).unwrap(), b"first");
        write_atomic(&output, b"second", true).unwrap();
        assert_eq!(std::fs::read(&output).unwrap(), b"second");
        let _ = std::fs::remove_file(output);
    }

    #[test]
    fn dropping_uncommitted_file_cleans_the_temporary_path() {
        let output = path("drop");
        let temporary = {
            let mut file = AtomicFile::create(&output).unwrap();
            file.file_mut().write_all(b"partial").unwrap();
            file.temporary_path().to_path_buf()
        };
        assert!(!temporary.exists());
        assert!(!output.exists());
    }
}
