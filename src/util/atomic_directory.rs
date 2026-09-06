//! Create-new publication of a complete artifact directory on supported hosts.

use std::fs;
use std::io;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

static NEXT: AtomicU64 = AtomicU64::new(0);

/// A sibling staging directory. Children must be fully written and synced before
/// commit. Existing destinations, including empty directories, are never replaced.
#[derive(Debug)]
pub struct AtomicDirectory {
    staging: PathBuf,
    destination: PathBuf,
    committed: bool,
}

impl AtomicDirectory {
    /// Reserve a unique staging directory without creating the successful name.
    pub fn create(destination: &Path) -> io::Result<Self> {
        match fs::symlink_metadata(destination) {
            Ok(_) => {
                return Err(io::Error::new(
                    io::ErrorKind::AlreadyExists,
                    "export destination already exists",
                ))
            }
            Err(e) if e.kind() == io::ErrorKind::NotFound => (),
            Err(e) => return Err(e),
        }
        let parent = destination
            .parent()
            .filter(|p| !p.as_os_str().is_empty())
            .unwrap_or_else(|| Path::new("."));
        for _ in 0..128 {
            let staging = parent.join(format!(
                ".rosalind-dataset-{}-{}.partial",
                std::process::id(),
                NEXT.fetch_add(1, Ordering::Relaxed)
            ));
            match fs::create_dir(&staging) {
                Ok(()) => {
                    return Ok(Self {
                        staging,
                        destination: destination.to_owned(),
                        committed: false,
                    })
                }
                Err(e) if e.kind() == io::ErrorKind::AlreadyExists => continue,
                Err(e) => return Err(e),
            }
        }
        Err(io::Error::new(
            io::ErrorKind::AlreadyExists,
            "cannot reserve export staging directory",
        ))
    }

    /// Path for incomplete private files, removed automatically on abandonment.
    pub fn path(&self) -> &Path {
        &self.staging
    }

    /// Atomically make the complete directory visible without replacing a raced
    /// destination. Publication requires Linux or macOS's exclusive rename API.
    pub fn commit(mut self) -> io::Result<PathBuf> {
        fs::File::open(&self.staging)?.sync_all()?;
        rename_exclusive(&self.staging, &self.destination)?;
        self.committed = true;
        if let Some(parent) = self.destination.parent() {
            if let Ok(directory) = fs::File::open(parent) {
                let _ = directory.sync_all();
            }
        }
        Ok(self.destination.clone())
    }
}

impl Drop for AtomicDirectory {
    fn drop(&mut self) {
        if !self.committed {
            let _ = fs::remove_dir_all(&self.staging);
        }
    }
}

#[cfg(any(target_os = "linux", target_os = "macos"))]
fn rename_exclusive(source: &Path, target: &Path) -> io::Result<()> {
    use std::ffi::CString;
    use std::os::unix::ffi::OsStrExt;
    let source = CString::new(source.as_os_str().as_bytes())
        .map_err(|_| io::Error::new(io::ErrorKind::InvalidInput, "NUL in staging path"))?;
    let target = CString::new(target.as_os_str().as_bytes())
        .map_err(|_| io::Error::new(io::ErrorKind::InvalidInput, "NUL in destination path"))?;
    // SAFETY: both NUL-terminated paths remain live through the synchronous call.
    // The no-replace flag prevents a check/rename race from clobbering a directory.
    #[cfg(target_os = "linux")]
    let result = unsafe {
        libc::syscall(
            libc::SYS_renameat2,
            libc::AT_FDCWD,
            source.as_ptr(),
            libc::AT_FDCWD,
            target.as_ptr(),
            libc::RENAME_NOREPLACE,
        )
    };
    #[cfg(target_os = "macos")]
    let result = unsafe { libc::renamex_np(source.as_ptr(), target.as_ptr(), libc::RENAME_EXCL) };
    if result == 0 {
        Ok(())
    } else {
        Err(io::Error::last_os_error())
    }
}

#[cfg(not(any(target_os = "linux", target_os = "macos")))]
fn rename_exclusive(_: &Path, _: &Path) -> io::Result<()> {
    Err(io::Error::new(
        io::ErrorKind::Unsupported,
        "atomic dataset directory publication requires Linux or macOS",
    ))
}

#[cfg(all(test, any(target_os = "linux", target_os = "macos")))]
mod tests {
    use super::*;

    #[test]
    fn publish_is_complete_and_raced_empty_destination_is_preserved() {
        let root = std::env::temp_dir().join(format!(
            "rosalind-directory-test-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&root).unwrap();
        let target = root.join("export");
        let stage = AtomicDirectory::create(&target).unwrap();
        fs::write(stage.path().join("part.parquet"), b"complete").unwrap();
        assert!(!target.exists());
        stage.commit().unwrap();
        assert_eq!(fs::read(target.join("part.parquet")).unwrap(), b"complete");
        assert!(AtomicDirectory::create(&target).is_err());
        let raced = root.join("raced");
        let stage = AtomicDirectory::create(&raced).unwrap();
        let abandoned = stage.path().to_owned();
        fs::write(stage.path().join("partial"), b"do not publish").unwrap();
        fs::create_dir(&raced).unwrap();
        assert!(stage.commit().is_err());
        assert!(raced.is_dir());
        assert_eq!(fs::read_dir(&raced).unwrap().count(), 0);
        assert!(!abandoned.exists());
        fs::remove_dir_all(root).unwrap();
    }
}
