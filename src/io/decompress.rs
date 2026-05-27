//! Transparent decompression + stream plumbing for pipe-native IO.
//!
//! Input may be plain text or gzip-compressed; bgzf (the BAM/"block gzip"
//! format) is a series of concatenated gzip members, so a multi-member gzip
//! decoder reads it transparently for sequential access. Random-access bgzf
//! (virtual offsets, `.gzi`/`.csi`) is Phase D and is not handled here.

use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Write};
use std::path::Path;

use flate2::bufread::MultiGzDecoder;

/// The two-byte gzip magic prefix (`1f 8b`), shared by plain gzip and bgzf.
const GZIP_MAGIC: [u8; 2] = [0x1f, 0x8b];

/// Wrap a buffered reader, transparently decompressing if it begins with the
/// gzip magic (`1f 8b`). Plain input passes through unchanged.
///
/// The two magic bytes are read out and then restored via [`Read::chain`], so
/// the returned reader still yields the complete original stream. The probe
/// loops because [`Read::read`] (like `fill_buf`) may return fewer bytes than
/// asked for — this keeps detection correct even for sources that deliver one
/// byte at a time (e.g. a trickling pipe).
pub fn maybe_decompress<R: BufRead + 'static>(mut reader: R) -> io::Result<Box<dyn BufRead>> {
    let mut magic = [0u8; 2];
    let mut filled = 0;
    while filled < magic.len() {
        match reader.read(&mut magic[filled..]) {
            Ok(0) => break, // EOF before two bytes: cannot be gzip
            Ok(n) => filled += n,
            Err(e) => return Err(e),
        }
    }
    let restored = io::Cursor::new(magic[..filled].to_vec()).chain(reader);
    if filled == 2 && magic == GZIP_MAGIC {
        // `MultiGzDecoder` decodes concatenated members (covers bgzf) but is
        // `Read`-only, so wrap it back into a `BufReader` for the `BufRead` return.
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(restored))))
    } else {
        Ok(Box::new(restored))
    }
}

/// Open `path` (or `-` for stdin) as a transparently-decompressed buffered reader.
pub fn open_input<P: AsRef<Path>>(path: P) -> io::Result<Box<dyn BufRead>> {
    let path = path.as_ref();
    if path.as_os_str() == "-" {
        maybe_decompress(BufReader::new(io::stdin()))
    } else {
        maybe_decompress(BufReader::new(File::open(path)?))
    }
}

/// Open `path` (or `-` for stdout) as a buffered writer.
pub fn open_output<P: AsRef<Path>>(path: P) -> io::Result<Box<dyn Write>> {
    let path = path.as_ref();
    if path.as_os_str() == "-" {
        Ok(Box::new(io::BufWriter::new(io::stdout())))
    } else {
        Ok(Box::new(io::BufWriter::new(File::create(path)?)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::write::GzEncoder;
    use flate2::Compression;
    use std::io::Cursor;

    fn gzip(bytes: &[u8]) -> Vec<u8> {
        let mut enc = GzEncoder::new(Vec::new(), Compression::default());
        enc.write_all(bytes).unwrap();
        enc.finish().unwrap()
    }

    #[test]
    fn plain_input_passes_through() {
        let mut out = String::new();
        maybe_decompress(Cursor::new(b"plain text".to_vec()))
            .unwrap()
            .read_to_string(&mut out)
            .unwrap();
        assert_eq!(out, "plain text");
    }

    #[test]
    fn gzip_input_is_transparently_decompressed() {
        let gz = gzip(b"hello world");
        let mut out = String::new();
        maybe_decompress(Cursor::new(gz))
            .unwrap()
            .read_to_string(&mut out)
            .unwrap();
        assert_eq!(out, "hello world");
    }

    #[test]
    fn empty_input_is_treated_as_plain_and_yields_nothing() {
        let mut out = String::new();
        maybe_decompress(Cursor::new(Vec::<u8>::new()))
            .unwrap()
            .read_to_string(&mut out)
            .unwrap();
        assert_eq!(out, "");
    }

    /// A `BufRead` that releases at most one byte per call, to prove magic
    /// detection survives sources that don't deliver two bytes at once.
    struct Trickle {
        data: Vec<u8>,
        pos: usize,
    }

    impl Read for Trickle {
        fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
            if self.pos >= self.data.len() || buf.is_empty() {
                return Ok(0);
            }
            buf[0] = self.data[self.pos];
            self.pos += 1;
            Ok(1)
        }
    }

    impl BufRead for Trickle {
        fn fill_buf(&mut self) -> io::Result<&[u8]> {
            let end = (self.pos + 1).min(self.data.len());
            Ok(&self.data[self.pos..end])
        }
        fn consume(&mut self, amt: usize) {
            self.pos += amt;
        }
    }

    #[test]
    fn gzip_detected_even_when_source_yields_one_byte_at_a_time() {
        let gz = gzip(b"trickled gzip payload");
        let reader = Trickle { data: gz, pos: 0 };
        let mut out = String::new();
        maybe_decompress(reader)
            .unwrap()
            .read_to_string(&mut out)
            .unwrap();
        assert_eq!(out, "trickled gzip payload");
    }
}
