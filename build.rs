fn main() {
    // These files are compiled into the Receipt Studio binary. Explicit change
    // tracking prevents a restored Cargo cache from serving stale verifier
    // assets when only the generated web bundle changed.
    println!("cargo:rerun-if-changed=web/verify/index.html");
    println!("cargo:rerun-if-changed=web/verify/pkg/rosalind_verify.js");
    println!("cargo:rerun-if-changed=web/verify/pkg/rosalind_verify_bg.wasm");
    rosalind_build_info::emit();
}
