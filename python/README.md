# Python Bindings for Rosalind

## Installation

### Prerequisites
- Python 3.8+
- Rust toolchain (install via [rustup](https://rustup.rs))
- maturin: `pip install maturin`

### Build and Install
```bash
cd /path/to/rosalind
maturin develop
# Or for a release build:
maturin develop --release
```

## Usage
```python
from rosalind_py import PyGenomicEngine

engine = PyGenomicEngine()
print(engine.list_plugins())

# region_start, region_end, reads[(pos, seq)], block_size
coverage = engine.run_rna_seq_plugin(
    region_start=100_000,
    region_end=101_000,
    reads=[(100_020, "ACGTACGT"), (100_050, "TTTACGT")],
    block_size=512,
)
print(coverage[:5])
```

Refer to the main project [README](../README.md) for more details about the available plugins and workloads.
