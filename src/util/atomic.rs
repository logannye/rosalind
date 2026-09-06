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

    /// Close the reservation handle while keeping the temporary path. This lets
    /// path-based writers such as htslib create/truncate the reserved sibling file.
    pub fn close_for_path_writer(&mut self) {
        self.file.take();
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
        } else {
            OpenOptions::new()
                .write(true)
                .open(&self.temporary)?
                .sync_all()?;
        }
        if replace {
            #[cfg(windows)]
            if destination.exists() {
                std::fs::remove_file(destination)?;
            }
            std::fs::rename(&self.temporary, destination)?;
        } else {
            // A same-filesystem hard link is an atomic create-new publication: it
            // cannot overwrite a destination that appears between preflight and
            // commit. Unlinking the temporary name leaves the committed inode.
            std::fs::hard_link(&self.temporary, destination)?;
            std::fs::remove_file(&self.temporary)?;
        }
        self.committed = true;
        Ok(destination.to_path_buf())
    }
}

/// Refuse an existing final destination before computation unless replacement was
/// explicitly requested.
pub fn ensure_destination(path: &Path, replace: bool) -> io::Result<()> {
    if !replace && path.exists() {
        Err(io::Error::new(
            io::ErrorKind::AlreadyExists,
            format!("destination already exists: {}", path.display()),
        ))
    } else {
        Ok(())
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

/// Publish a small group of staged artifacts, rolling back completed publications
/// if a later destination fails. Existing files are preserved with same-directory
/// hard links before replacement. This handles reported I/O failures; it does not
/// promise a multi-file atomic view to concurrent readers or across process death.
pub fn commit_group(files: Vec<(AtomicFile, PathBuf)>, replace: bool) -> io::Result<()> {
    let mut backups = Vec::with_capacity(files.len());
    for (_, destination) in &files {
        let backup = if replace && destination.exists() {
            let mut backup = AtomicFile::create(destination)?;
            backup.close_for_path_writer();
            std::fs::remove_file(backup.temporary_path())?;
            std::fs::hard_link(destination, backup.temporary_path())?;
            Some(backup)
        } else {
            None
        };
        backups.push(backup);
    }
    let mut published: Vec<(usize, PathBuf)> = Vec::with_capacity(files.len());
    for (index, (file, destination)) in files.into_iter().enumerate() {
        if let Err(error) = file.commit_as(&destination, replace) {
            let mut rollback_error = None;
            for (prior_index, prior_path) in published.iter().rev() {
                let result = if let Some(backup) = &backups[*prior_index] {
                    std::fs::rename(backup.temporary_path(), prior_path)
                } else {
                    std::fs::remove_file(prior_path)
                };
                if let Err(error) = result {
                    rollback_error = Some(error);
                }
            }
            return Err(match rollback_error {
                Some(rollback) => io::Error::other(format!(
                    "publication failed: {error}; rollback also failed: {rollback}"
                )),
                None => error,
            });
        }
        published.push((index, destination));
    }
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

    #[test]
    fn group_failure_restores_replaced_outputs_and_leaves_no_new_success() {
        let first = path("group-existing");
        let fresh = path("group-fresh");
        let last = path("group-invalid");
        write_atomic(&first, b"original", false).unwrap();
        let stage = |path: &Path, bytes: &[u8]| {
            let mut file = AtomicFile::create(path).unwrap();
            file.file_mut().write_all(bytes).unwrap();
            file
        };
        let files = vec![
            (stage(&first, b"replacement"), first.clone()),
            (stage(&fresh, b"new"), fresh.clone()),
            (stage(&last, b"last"), last.clone()),
        ];
        // A directory cannot be replaced by the staged regular file. It appears
        // after backup preparation would normally occur; a missing parent is an
        // equally useful failure injected at the last commit.
        let invalid = path("group-missing-parent").join("last");
        let mut files = files;
        files[2].1 = invalid;
        assert!(commit_group(files, true).is_err());
        assert_eq!(std::fs::read(&first).unwrap(), b"original");
        assert!(!fresh.exists());
        let _ = std::fs::remove_file(first);
    }
}
