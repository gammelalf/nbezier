use criterion::{black_box, criterion_group, criterion_main, Criterion};
mod custom_profiler;
mod samples;
use samples::CURVES;

pub fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("split", |b| {
        for curve in CURVES.iter() {
            b.iter(|| black_box(curve.split(0.5)))
        }
    });
    c.bench_function("eval", |b| {
        for curve in CURVES.iter() {
            b.iter(|| black_box(curve.castlejau_eval(0.5)))
        }
    });
}

criterion_group!{
    name = benches;
    config = Criterion::default().with_profiler(custom_profiler::FlamegraphProfiler::new(100));
    targets = criterion_benchmark
}
criterion_main!(benches);