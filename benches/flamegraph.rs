use criterion::{criterion_group, criterion_main, Criterion};
use pprof::criterion::{PProfProfiler, Output};
mod common;

criterion_group!{
    name = benches;
    config = Criterion::default().with_profiler(PProfProfiler::new(100, Output::Flamegraph(None)));
    targets = common::bezier::all
}
criterion_main!(benches);