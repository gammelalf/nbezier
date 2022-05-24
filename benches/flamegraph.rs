use criterion::{criterion_group, criterion_main, Criterion};
mod common;
use common::flamegraph_profiler::FlamegraphProfiler;

criterion_group!{
    name = benches;
    config = Criterion::default().with_profiler(FlamegraphProfiler::new(100));
    targets = common::bezier::all
}
criterion_main!(benches);