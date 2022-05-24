use criterion::{black_box, Criterion};
use crate::common::samples::CURVES;

pub fn split(c: &mut Criterion) {
    c.bench_function("split", |b| {
        for curve in CURVES.iter() {
            b.iter(|| black_box(curve.split(0.5)))
        }
    });
}

pub fn eval(c: &mut Criterion) {
    c.bench_function("eval", |b| {
        for curve in CURVES.iter() {
            b.iter(|| black_box(curve.castlejau_eval(0.5)))
        }
    });
}

pub fn normal(c: &mut Criterion) {
    c.bench_function("normal", |b| {
        for curve in CURVES.iter() {
            b.iter(|| black_box(curve.normal(0.5)))
        }
    });
}

pub fn all(c: &mut Criterion) {
    //split(c);
    //eval(c);
    normal(c);
}
