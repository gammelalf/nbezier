use std::time::Duration;
use criterion::{criterion_group, criterion_main, Criterion};
use pprof::criterion::{PProfProfiler, Output};
mod common;

static SECOND: Duration = Duration::from_secs(1);
criterion_group!{
    name = benches;
    config = Criterion::default()
    .with_profiler(PProfProfiler::new(100, Output::Flamegraph(None)))
    ;
    targets = common::bezier::all
}
criterion_main!(benches);