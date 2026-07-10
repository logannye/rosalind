use std::io::{BufRead, BufReader, Read, Write};
use std::net::TcpStream;
use std::process::{Command, Stdio};

use rosalind::provenance::RunManifest;

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_rosalind")
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

    let mut child = Command::new(bin())
        .args([
            "studio",
            receipt_path.to_str().unwrap(),
            "--no-open",
            "--port",
            "0",
            "--json",
        ])
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    let mut line = String::new();
    BufReader::new(child.stdout.take().unwrap())
        .read_line(&mut line)
        .unwrap();
    assert!(line.contains("\"bind\":\"127.0.0.1\""), "{line}");
    assert!(line.contains("\"preloaded\":1"), "{line}");
    let marker = "\"port\":";
    let start = line.find(marker).unwrap() + marker.len();
    let port = line[start..]
        .split(|character: char| !character.is_ascii_digit())
        .next()
        .unwrap()
        .parse::<u16>()
        .unwrap();

    let response = get(port, "/");
    assert!(response.starts_with("HTTP/1.1 200 OK"), "{response}");
    assert!(response.contains("Content-Security-Policy:"), "{response}");
    assert!(response.contains("rosalind-studio-test-"), "{response}");
    assert!(response.contains("Receipt Studio"), "{response}");

    let js = get(port, "/pkg/rosalind_verify.js");
    assert!(js.starts_with("HTTP/1.1 200 OK"), "{js}");
    assert!(
        js.contains("evaluate_trust"),
        "generated bindings are embedded"
    );

    child.kill().ok();
    child.wait().ok();
    std::fs::remove_file(receipt_path).ok();
}

fn get(port: u16, path: &str) -> String {
    let mut stream = TcpStream::connect(("127.0.0.1", port)).unwrap();
    write!(
        stream,
        "GET {path} HTTP/1.1\r\nHost: 127.0.0.1\r\nConnection: close\r\n\r\n"
    )
    .unwrap();
    let mut response = String::new();
    stream.read_to_string(&mut response).unwrap();
    response
}
