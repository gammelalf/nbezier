///
/// Use the pprof-rs profiler outputting a svg flamegraph.
///
/// Needs to be called with `--profile-time=<seconds>`
///
use criterion::{criterion_group, criterion_main, Criterion};
use pprof::criterion::{Output, PProfProfiler};
mod common;

criterion_group! {
    name = benches;
    config = Criterion::default()
    .with_profiler(PProfProfiler::new(100, Output::Flamegraph(None)))
    ;
    targets = common::bezier::all
}
criterion_main!(benches);
