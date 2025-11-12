//! Performance benchmarks

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use rosalind::*;

fn benchmark_simulation(c: &mut Criterion) {
    // TODO: Setup benchmark machine

    c.bench_function("simulate_t=10000", |b| {
        b.iter(|| {
            // TODO: Run simulation
            black_box(());
        });
    });
}

criterion_group!(benches, benchmark_simulation);
criterion_main!(benches);
