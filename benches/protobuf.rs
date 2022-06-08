///
/// Use the pprof-rs profiler outputting protobuf data.
///
/// Then use `go tool pprof` to render graphs from this data.
/// For example:
/// - `go tool pprof -svg profile.pb` for a flowchart
/// - `go tool pprof -http 127.0.0.1:8000` for an interactive view including flowchart and flamegraph
///
/// Needs to be called with `--profile-time=<seconds>`
///

use criterion::{criterion_group, criterion_main, Criterion};
use pprof::criterion::{PProfProfiler, Output};
mod common;

criterion_group!{
    name = benches;
    config = Criterion::default()
    .with_profiler(PProfProfiler::new(100, Output::Protobuf))
    ;
    targets = common::bezier::all
}
criterion_main!(benches);