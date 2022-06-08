use std::time::Duration;
use criterion::{criterion_group, criterion_main, Criterion};
mod common;

static SECOND: Duration = Duration::from_secs(1);
criterion_group!{
    name = benches;
    config = Criterion::default()
    .warm_up_time(SECOND)
    .measurement_time(SECOND)
    ;
    targets = common::bezier::all
}
criterion_main!(benches);
