//! Explicit pair-table input, bounded and guarded inside the artifact lifetime.
use super::descriptor::CohortLimits;
use super::pairs::{PairLimits, PairSpec};
use super::{CohortError, Result};
use crate::dataset::InputSnapshot;
use std::fs::{self, File};
use std::io::Read;
use std::path::PathBuf;

#[derive(Debug, Clone)]
pub(crate) struct PairInput {
    pub path: PathBuf,
    pub max_table_bytes: usize,
    pub limits: PairLimits,
}

pub(crate) struct LoadedPairs {
    pub path: PathBuf,
    pub blake3: String,
    pub guard: InputSnapshot,
    pub pairs: Vec<PairSpec>,
    pub limits: PairLimits,
}

pub(crate) fn read_pairs(input: &PairInput, limits: CohortLimits) -> Result<LoadedPairs> {
    if input.max_table_bytes == 0 || input.limits.max_pairs == 0 {
        return Err(CohortError::Limit(
            "pair table envelopes must be positive".into(),
        ));
    }
    let path = fs::canonicalize(&input.path)?;
    let guard = InputSnapshot::capture([path.clone()])?;
    let length = usize::try_from(fs::metadata(&path)?.len())
        .map_err(|_| CohortError::Limit("pair table exceeds address space".into()))?;
    if length > input.max_table_bytes {
        return Err(CohortError::Limit(
            "pair table exceeds its byte envelope".into(),
        ));
    }
    // Text, parsed records and scope resolution coexist briefly. Admit before
    // allocating, regardless of whether this is an execution or a metadata plan.
    limits.admit((length as u64).saturating_mul(32).saturating_add(128 << 10))?;
    let mut bytes = Vec::with_capacity(length);
    File::open(&path)?
        .take((length as u64).saturating_add(1))
        .read_to_end(&mut bytes)?;
    guard.verify()?;
    if bytes.len() != length {
        return Err(CohortError::Corrupt(
            "pair table changed during reading".into(),
        ));
    }
    let hash = blake3::hash(&bytes).to_hex().to_string();
    let text = std::str::from_utf8(&bytes)
        .map_err(|_| CohortError::Incompatible("pair table must be UTF-8 TSV".into()))?;
    let mut lines = text.lines();
    if lines.next() != Some("id\tleft\tright") {
        return Err(CohortError::Incompatible(
            "pair table header must be id<TAB>left<TAB>right".into(),
        ));
    }
    let mut pairs = Vec::new();
    let mut metadata = 4096u64;
    for line in lines {
        limits.admit(0)?;
        if pairs.len() >= input.limits.max_pairs {
            return Err(CohortError::Limit(
                "pair table exceeds its pair-count envelope".into(),
            ));
        }
        // Validate bounded text before owned strings and subsequent lookup state.
        let mut columns = line.split('\t');
        let mut next = || -> Result<String> {
            let value = columns.next().ok_or_else(|| {
                CohortError::Incompatible("pair table requires exactly three columns".into())
            })?;
            if value.is_empty() || value.len() > 256 || value.chars().any(char::is_control) {
                return Err(CohortError::Incompatible("pair IDs and member IDs must be nonempty, at most 256 UTF-8 bytes, without controls".into()));
            }
            Ok(value.to_owned())
        };
        let pair = PairSpec {
            id: next()?,
            left: next()?,
            right: next()?,
        };
        if columns.next().is_some() {
            return Err(CohortError::Incompatible(
                "pair table requires exactly three columns".into(),
            ));
        }
        metadata = metadata
            .saturating_add(1024)
            .saturating_add((pair.id.len() + pair.left.len() + pair.right.len()) as u64 * 4);
        if metadata > input.limits.max_metadata_bytes {
            return Err(CohortError::Limit(
                "pair metadata exceeds its byte envelope".into(),
            ));
        }
        limits.admit(metadata)?;
        pairs.push(pair);
    }
    guard.verify()?;
    Ok(LoadedPairs {
        path,
        blake3: hash,
        guard,
        pairs,
        limits: input.limits,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn table_parser_bounds_cardinality_text_and_preserves_explicit_order() {
        let root = std::env::temp_dir().join(format!("rosalind-pair-table-{}", std::process::id()));
        fs::create_dir(&root).unwrap();
        struct Cleanup(PathBuf);
        impl Drop for Cleanup {
            fn drop(&mut self) {
                let _ = fs::remove_dir_all(&self.0);
            }
        }
        let _cleanup = Cleanup(root.clone());
        let path = root.join("pairs.tsv");
        let input = PairInput {
            path: path.clone(),
            max_table_bytes: 1024,
            limits: PairLimits {
                max_pairs: 2,
                ..PairLimits::default()
            },
        };
        fs::write(&path, "id\tleft\tright\nz\tB\tA\na\tA\tB\n").unwrap();
        let loaded = read_pairs(&input, CohortLimits::default()).unwrap();
        assert_eq!(loaded.pairs[0].id, "z");
        assert_eq!(loaded.pairs[0].left, "B");
        assert_eq!(loaded.pairs[1].id, "a");
        assert_eq!(
            loaded.blake3,
            blake3::hash(&fs::read(&path).unwrap()).to_hex().to_string()
        );
        fs::write(&path, "id\tleft\tright\nz\tB\tA\na\tA\tB\nextra\tA\tB\n").unwrap();
        assert!(matches!(
            read_pairs(&input, CohortLimits::default()),
            Err(CohortError::Limit(_))
        ));
        assert!(loaded.guard.verify().is_err());
        for invalid in [
            b"id\tleft\tright\n\xff\tA\tB\n".as_slice(),
            b"id\tleft\tright\nx\tA\tB\textra\n",
            b"id\tleft\tright\n\n",
            b"id\tleft\tright\nx\t\tB\n",
        ] {
            fs::write(&path, invalid).unwrap();
            assert!(read_pairs(&input, CohortLimits::default()).is_err());
        }
        fs::write(&path, "id\tleft\tright\n").unwrap();
        assert!(read_pairs(&input, CohortLimits::default())
            .unwrap()
            .pairs
            .is_empty());
        let tiny = PairInput {
            max_table_bytes: 1,
            ..input
        };
        assert!(matches!(
            read_pairs(&tiny, CohortLimits::default()),
            Err(CohortError::Limit(_))
        ));
    }
}
