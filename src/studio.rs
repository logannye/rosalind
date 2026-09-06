//! Loopback-only server for the embedded, client-side Receipt Studio.

use std::io::{Read, Write};
use std::net::{TcpListener, TcpStream};
use std::path::PathBuf;
use std::process::Command;
use std::time::{Duration, Instant};

use anyhow::{Context, Result};

const INDEX_HTML: &str = include_str!("../web/verify/index.html");
const STUDIO_JS: &[u8] = include_bytes!("../web/verify/pkg/rosalind_verify.js");
const STUDIO_WASM: &[u8] = include_bytes!("../web/verify/pkg/rosalind_verify_bg.wasm");
const MAX_REQUEST_HEADER_BYTES: usize = 8192;
const REQUEST_HEADER_TIMEOUT: Duration = Duration::from_secs(5);

/// Local Studio launch configuration.
#[derive(Debug, Clone)]
pub struct StudioSpec {
    /// Receipts and certificates preloaded into the page.
    pub receipts: Vec<PathBuf>,
    /// Suppress opening the system browser.
    pub no_open: bool,
    /// Loopback port, or zero to ask the OS for an unused port.
    pub port: u16,
    /// Print one startup JSON object rather than a human message.
    pub json: bool,
}

/// Serve embedded Studio assets on `127.0.0.1` until interrupted.
pub fn serve_studio(spec: &StudioSpec) -> Result<()> {
    let listener = TcpListener::bind(("127.0.0.1", spec.port))
        .with_context(|| format!("failed to bind Receipt Studio on port {}", spec.port))?;
    let address = listener.local_addr()?;
    let url = format!("http://127.0.0.1:{}/", address.port());
    let html = inject_preloads(INDEX_HTML, &spec.receipts)?;
    if spec.json {
        println!(
            "{{\"schema\":1,\"url\":\"{}\",\"bind\":\"127.0.0.1\",\"port\":{},\"preloaded\":{}}}",
            url,
            address.port(),
            spec.receipts.len()
        );
    } else {
        println!("Receipt Studio: {url}");
        println!("Browser-local artifacts only; press Ctrl-C to stop.");
    }
    std::io::stdout().flush()?;
    if !spec.no_open {
        open_browser(&url);
    }
    for stream in listener.incoming() {
        match stream {
            Ok(mut stream) => {
                if let Err(error) = respond(&mut stream, html.as_bytes()) {
                    eprintln!("studio request failed: {error}");
                }
            }
            Err(error) => eprintln!("studio connection failed: {error}"),
        }
    }
    Ok(())
}

fn inject_preloads(template: &str, paths: &[PathBuf]) -> Result<String> {
    const EMPTY: &str = "<script type=\"application/json\" id=\"studio-preload\">[]</script>";
    let mut entries = Vec::with_capacity(paths.len());
    for path in paths {
        let text = std::fs::read_to_string(path)
            .with_context(|| format!("failed to preload receipt {}", path.display()))?;
        let name = path
            .file_name()
            .and_then(|name| name.to_str())
            .unwrap_or("receipt.json");
        entries.push(format!(
            "{{\"name\":\"{}\",\"text\":\"{}\"}}",
            json_for_script(name),
            json_for_script(&text)
        ));
    }
    let replacement = format!(
        "<script type=\"application/json\" id=\"studio-preload\">[{}]</script>",
        entries.join(",")
    );
    Ok(template.replacen(EMPTY, &replacement, 1))
}

fn respond(stream: &mut TcpStream, html: &[u8]) -> std::io::Result<()> {
    let started = Instant::now();
    let mut request = [0_u8; MAX_REQUEST_HEADER_BYTES];
    let length = read_request_headers(&mut request, |buffer| {
        // Bound the whole header exchange, including clients that trickle a
        // byte at a time. A fresh full timeout per read would not do that.
        let remaining = REQUEST_HEADER_TIMEOUT
            .checked_sub(started.elapsed())
            .filter(|remaining| !remaining.is_zero())
            .ok_or_else(|| {
                std::io::Error::new(std::io::ErrorKind::TimedOut, "request header timed out")
            })?;
        stream.set_read_timeout(Some(remaining))?;
        stream.read(buffer)
    })?;
    stream.set_write_timeout(Some(REQUEST_HEADER_TIMEOUT))?;
    respond_to_request(stream, &request[..length], html)
}

/// Assemble headers within the caller's fixed buffer. TCP read boundaries do
/// not delimit requests, even when a client makes a single formatted write.
fn read_request_headers(
    request: &mut [u8],
    mut read: impl FnMut(&mut [u8]) -> std::io::Result<usize>,
) -> std::io::Result<usize> {
    let mut used: usize = 0;
    loop {
        if used == request.len() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "request header exceeds 8192 bytes",
            ));
        }
        let count = match read(&mut request[used..]) {
            Err(error) if error.kind() == std::io::ErrorKind::Interrupted => continue,
            result => result?,
        };
        if count == 0 {
            return Err(std::io::Error::new(
                std::io::ErrorKind::UnexpectedEof,
                "connection closed before complete request headers",
            ));
        }
        let search_start = used.saturating_sub(3);
        used += count;
        if let Some(end) = request[search_start..used]
            .windows(4)
            .position(|bytes| bytes == b"\r\n\r\n")
        {
            return Ok(search_start + end + 4);
        }
    }
}

fn respond_to_request(stream: &mut impl Write, request: &[u8], html: &[u8]) -> std::io::Result<()> {
    let malformed = || {
        std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "malformed HTTP request line",
        )
    };
    let line_end = request
        .windows(2)
        .position(|bytes| bytes == b"\r\n")
        .ok_or_else(malformed)?;
    let line = std::str::from_utf8(&request[..line_end]).map_err(|_| malformed())?;
    let mut parts = line.split_ascii_whitespace();
    let _method = parts.next().ok_or_else(malformed)?;
    let target = parts.next().ok_or_else(malformed)?;
    let version = parts.next().ok_or_else(malformed)?;
    if parts.next().is_some() || !matches!(version, "HTTP/1.0" | "HTTP/1.1") {
        return Err(malformed());
    }
    let path = target.split('?').next().ok_or_else(malformed)?;
    let (status, content_type, body): (&str, &str, &[u8]) = match path {
        "/" | "/index.html" => ("200 OK", "text/html; charset=utf-8", html),
        "/pkg/rosalind_verify.js" => ("200 OK", "text/javascript; charset=utf-8", STUDIO_JS),
        "/pkg/rosalind_verify_bg.wasm" => ("200 OK", "application/wasm", STUDIO_WASM),
        _ => ("404 Not Found", "text/plain; charset=utf-8", b"not found\n"),
    };
    write!(
        stream,
        "HTTP/1.1 {status}\r\nContent-Type: {content_type}\r\nContent-Length: {}\r\nCache-Control: no-store\r\nContent-Security-Policy: default-src 'self'; script-src 'self' 'unsafe-inline' 'wasm-unsafe-eval'; style-src 'self' 'unsafe-inline'; img-src 'self' data:; connect-src 'self'\r\nX-Content-Type-Options: nosniff\r\nConnection: close\r\n\r\n",
        body.len()
    )?;
    stream.write_all(body)
}

fn open_browser(url: &str) {
    #[cfg(target_os = "macos")]
    let mut command = Command::new("open");
    #[cfg(all(unix, not(target_os = "macos")))]
    let mut command = Command::new("xdg-open");
    #[cfg(windows)]
    let mut command = {
        let mut command = Command::new("cmd");
        command.args(["/C", "start", ""]);
        command
    };
    let _ = command.arg(url).spawn();
}

fn json_for_script(value: &str) -> String {
    let mut output = String::with_capacity(value.len());
    for character in value.chars() {
        match character {
            '"' => output.push_str("\\\""),
            '\\' => output.push_str("\\\\"),
            '\n' => output.push_str("\\n"),
            '\r' => output.push_str("\\r"),
            '\t' => output.push_str("\\t"),
            '<' => output.push_str("\\u003c"),
            '>' => output.push_str("\\u003e"),
            '&' => output.push_str("\\u0026"),
            value if (value as u32) < 0x20 => {
                output.push_str(&format!("\\u{:04x}", value as u32));
            }
            value => output.push(value),
        }
    }
    output
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn one_byte_reads_assemble_the_complete_javascript_request() {
        let request = b"GET /pkg/rosalind_verify.js?cache=0 HTTP/1.1\r\nHost: localhost\r\n\r\n";
        let mut remaining = request.as_slice();
        let mut assembled = [0_u8; MAX_REQUEST_HEADER_BYTES];
        let length =
            read_request_headers(&mut assembled, |buffer| remaining.read(&mut buffer[..1]))
                .unwrap();
        assert_eq!(&assembled[..length], request);
        assert!(remaining.is_empty());

        let mut response = Vec::new();
        respond_to_request(&mut response, &assembled[..length], b"HTML fallback").unwrap();
        let split = response
            .windows(4)
            .position(|bytes| bytes == b"\r\n\r\n")
            .unwrap();
        let headers = std::str::from_utf8(&response[..split]).unwrap();
        assert!(headers.starts_with("HTTP/1.1 200 OK\r\n"));
        assert!(headers.contains("Content-Type: text/javascript; charset=utf-8"));
        assert!(headers.contains(&format!("Content-Length: {}\r\n", STUDIO_JS.len())));
        assert_eq!(&response[split + 4..], STUDIO_JS);
    }

    #[test]
    fn incomplete_and_oversized_headers_are_rejected_with_bounded_reads() {
        let mut assembled = [0_u8; MAX_REQUEST_HEADER_BYTES];
        let mut incomplete = b"GET /pkg/rosalind_verify.js HTTP/1.1\r\n".as_slice();
        let error =
            read_request_headers(&mut assembled, |buffer| incomplete.read(buffer)).unwrap_err();
        assert_eq!(error.kind(), std::io::ErrorKind::UnexpectedEof);

        let mut bytes_read = 0;
        let error = read_request_headers(&mut assembled, |buffer| {
            buffer.fill(b'x');
            bytes_read += buffer.len();
            Ok(buffer.len())
        })
        .unwrap_err();
        assert_eq!(error.kind(), std::io::ErrorKind::InvalidData);
        assert_eq!(bytes_read, MAX_REQUEST_HEADER_BYTES);

        let mut response = Vec::new();
        let error =
            respond_to_request(&mut response, b"GET \r\n\r\n", b"HTML fallback").unwrap_err();
        assert_eq!(error.kind(), std::io::ErrorKind::InvalidData);
        assert!(response.is_empty());
    }

    #[test]
    fn complete_headers_at_the_capacity_boundary_are_accepted() {
        let mut input = [b'x'; MAX_REQUEST_HEADER_BYTES];
        input[MAX_REQUEST_HEADER_BYTES - 4..].copy_from_slice(b"\r\n\r\n");
        let mut remaining = input.as_slice();
        let mut assembled = [0_u8; MAX_REQUEST_HEADER_BYTES];
        assert_eq!(
            read_request_headers(&mut assembled, |buffer| remaining.read(buffer)).unwrap(),
            MAX_REQUEST_HEADER_BYTES
        );
    }

    #[test]
    fn preload_json_cannot_close_its_script_element() {
        let path = std::env::temp_dir().join(format!(
            "rosalind-studio-preload-{}.json",
            std::process::id()
        ));
        std::fs::write(&path, "{\"path\":\"</script>\"}").unwrap();
        let html = inject_preloads(INDEX_HTML, std::slice::from_ref(&path)).unwrap();
        let preload = html
            .split("id=\"studio-preload\">")
            .nth(1)
            .unwrap()
            .split("</script>")
            .next()
            .unwrap();
        assert!(!preload.contains("</script>"));
        assert!(preload.contains("\\u003c/script\\u003e"));
        std::fs::remove_file(path).ok();
    }
}
