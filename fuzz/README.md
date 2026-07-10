# Receipt parser fuzz targets

The three cargo-fuzz targets exercise the canonical receipt parser, historical
human-readable command parser, and tokenized `command_argv` parser. They are kept
outside the release workspace and do not affect normal builds.

```sh
cargo install cargo-fuzz
cargo fuzz run canonical_json
cargo fuzz run legacy_command
cargo fuzz run command_argv
```
