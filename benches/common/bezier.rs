use criterion::{black_box, Criterion};
use gammalg::bezier::BezierCurve;
use crate::common::samples::CURVES;

macro_rules! bench_curve_types {
    ($group:ident, $body:expr) => {
        $group.bench_function("linear", |b| {
            for curve in CURVES.LINEAR.iter() {
                b.iter(|| $body(curve))
            }
        });
        $group.bench_function("quadratic", |b| {
            for curve in CURVES.QUADRATIC.iter() {
                b.iter(|| $body(curve))
            }
        });
        $group.bench_function("cubic", |b| {
            for curve in CURVES.CUBIC.iter() {
                b.iter(|| $body(curve))
            }
        });
        $group.bench_function("higher", |b| {
            for curve in CURVES.HIGHER.iter() {
                b.iter(|| $body(curve))
            }
        });
    }
}

pub fn split(c: &mut Criterion) {
    let mut g = c.benchmark_group("BezierCurve::split");
    bench_curve_types!(g, |curve: &BezierCurve<f64>| curve.split(black_box(0.5)));
}

pub fn castlejau_eval(c: &mut Criterion) {
    let mut g = c.benchmark_group("BezierCurve::castlejau_eval");
    bench_curve_types!(g, |curve: &BezierCurve<f64>| curve.castlejau_eval(black_box(0.5)));
}

pub fn normal(c: &mut Criterion) {
    let mut g = c.benchmark_group("BezierCurve::normal");
    bench_curve_types!(g, |curve: &BezierCurve<f64>| curve.normal(black_box(0.5)));
}

pub fn all(c: &mut Criterion) {
    split(c);
    castlejau_eval(c);
    normal(c);
}
