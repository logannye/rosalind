//! Loopback-only server for the embedded, client-side Receipt Studio.

use std::io::{Read, Write};
use std::net::{TcpListener, TcpStream};
use std::path::PathBuf;
use std::process::Command;

use anyhow::{Context, Result};

const INDEX_HTML: &str = include_str!("../web/verify/index.html");
const STUDIO_JS: &[u8] = include_bytes!("../web/verify/pkg/rosalind_verify.js");
const STUDIO_WASM: &[u8] = include_bytes!("../web/verify/pkg/rosalind_verify_bg.wasm");

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
    let mut request = [0_u8; 8192];
    let read = stream.read(&mut request)?;
    let line = String::from_utf8_lossy(&request[..read]);
    let path = line
        .lines()
        .next()
        .and_then(|line| line.split_whitespace().nth(1))
        .unwrap_or("/")
        .split('?')
        .next()
        .unwrap_or("/");
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
