use clap::Parser;

fn main() {
    let exit = xtask::run(xtask::Cli::parse());
    std::process::exit(exit);
}
