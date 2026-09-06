use std::io::{BufRead, BufReader, Read, Write};
use std::net::TcpStream;
use std::process::{Child, Command, Stdio};
use std::time::Duration;

use rosalind::provenance::RunManifest;

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_rosalind")
}

struct StudioChild(Child);

impl Drop for StudioChild {
    fn drop(&mut self) {
        self.0.kill().ok();
        self.0.wait().ok();
    }
}

fn start_studio(receipts: &[&std::path::Path]) -> (StudioChild, u16, String) {
    let mut child = StudioChild(
        Command::new(bin())
            .arg("studio")
            .args(receipts)
            .args(["--no-open", "--port", "0", "--json"])
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .unwrap(),
    );
    let mut line = String::new();
    BufReader::new(child.0.stdout.take().unwrap())
        .read_line(&mut line)
        .unwrap();
    assert!(line.contains("\"bind\":\"127.0.0.1\""), "{line}");
    let marker = "\"port\":";
    let start = line.find(marker).unwrap() + marker.len();
    let port = line[start..]
        .split(|character: char| !character.is_ascii_digit())
        .next()
        .unwrap()
        .parse::<u16>()
        .unwrap();
    (child, port, line)
}

#[test]
fn studio_binds_loopback_serves_embedded_assets_and_preloads_receipts() {
    let receipt_path = std::env::temp_dir().join(format!(
        "rosalind-studio-test-{}.manifest.json",
        std::process::id()
    ));
    let mut receipt = RunManifest::new("features");
    receipt.finalize();
    std::fs::write(&receipt_path, receipt.to_canonical_json()).unwrap();

    let (_child, port, line) = start_studio(&[&receipt_path]);
    assert!(line.contains("\"preloaded\":1"), "{line}");

    let response = get(port, "/");
    assert!(response.starts_with("HTTP/1.1 200 OK"), "{response}");
    assert!(response.contains("Content-Security-Policy:"), "{response}");
    assert!(response.contains("'wasm-unsafe-eval'"), "{response}");
    assert!(response.contains("rosalind-studio-test-"), "{response}");
    assert!(response.contains("Receipt Studio"), "{response}");
    assert!(response.contains("id=\"sample-receipts\""), "{response}");
    assert!(response.contains("Load sample chain"), "{response}");
    assert!(response.contains("not supplied"), "{response}");

    let js = get(port, "/pkg/rosalind_verify.js");
    assert!(js.starts_with("HTTP/1.1 200 OK"), "{js}");
    assert!(
        js.contains("evaluate_trust"),
        "generated bindings are embedded"
    );

    std::fs::remove_file(receipt_path).ok();
}

#[test]
fn studio_waits_for_fragmented_request_headers_before_selecting_the_asset() {
    let (_child, port, _) = start_studio(&[]);
    let mut stream = TcpStream::connect(("127.0.0.1", port)).unwrap();
    stream.set_nodelay(true).unwrap();
    stream
        .set_read_timeout(Some(Duration::from_millis(250)))
        .unwrap();

    // A formatted write can send "GET " separately from its path. The old
    // single-read implementation answered this fragment with the HTML page.
    stream.write_all(b"GET ").unwrap();
    assert_no_response_yet(&mut stream);
    stream
        .write_all(b"/pkg/rosalind_verify.js HTTP/1.1\r\nHost: 127.0.0.1\r\n")
        .unwrap();
    assert_no_response_yet(&mut stream);
    stream.write_all(b"Connection: close\r\n\r").unwrap();
    assert_no_response_yet(&mut stream);
    stream.write_all(b"\n").unwrap();
    stream
        .set_read_timeout(Some(Duration::from_secs(5)))
        .unwrap();

    let response = read_response(&mut stream);
    let (headers, body) = response.split_once("\r\n\r\n").unwrap();
    assert!(headers.starts_with("HTTP/1.1 200 OK"), "{headers}");
    assert!(headers.contains("Content-Type: text/javascript; charset=utf-8"));
    assert_eq!(
        body.as_bytes(),
        include_bytes!("../web/verify/pkg/rosalind_verify.js")
    );
}

fn assert_no_response_yet(stream: &mut TcpStream) {
    match stream.read(&mut [0_u8; 1]) {
        Err(error)
            if matches!(
                error.kind(),
                std::io::ErrorKind::WouldBlock | std::io::ErrorKind::TimedOut
            ) => {}
        result => panic!("server answered an incomplete request: {result:?}"),
    }
}

fn get(port: u16, path: &str) -> String {
    let mut stream = TcpStream::connect(("127.0.0.1", port)).unwrap();
    write!(
        stream,
        "GET {path} HTTP/1.1\r\nHost: 127.0.0.1\r\nConnection: close\r\n\r\n"
    )
    .unwrap();
    read_response(&mut stream)
}

fn read_response(stream: &mut TcpStream) -> String {
    let mut response = Vec::new();
    let mut chunk = [0u8; 16 * 1024];
    loop {
        match stream.read(&mut chunk) {
            Ok(0) => break,
            Ok(read) => response.extend_from_slice(&chunk[..read]),
            // macOS may report a reset when the tiny no-keepalive server closes
            // immediately after a complete response. The HTTP Content-Length
            // assertion below still proves whether the body was complete.
            Err(error) if error.kind() == std::io::ErrorKind::ConnectionReset => break,
            Err(error) => panic!("failed to read Studio response: {error}"),
        }
    }
    let response = String::from_utf8(response).unwrap();
    let (headers, body) = response.split_once("\r\n\r\n").unwrap();
    let length = headers
        .lines()
        .find_map(|line| line.strip_prefix("Content-Length: "))
        .unwrap()
        .parse::<usize>()
        .unwrap();
    assert_eq!(body.len(), length, "complete HTTP response body");
    response
}
