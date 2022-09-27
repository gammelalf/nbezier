use nbezier::nbezier::BezierCurve;
use nalgebra::{Matrix2x4, Matrix2xX};
use criterion::black_box;

#[global_allocator]
static ALLOC: dhat::Alloc = dhat::Alloc;

fn main() {
    let profiler = dhat::Profiler::new_heap();
    let curve = BezierCurve(Matrix2x4::new(
        50.0, 200.0, 0.0, 50.0, 0.0, 33.0, 66.0, 100.0
    ));
    black_box(curve.split(0.5));
    black_box(curve.castlejau_eval(0.5));
    black_box(curve);
    drop(profiler);

    let profiler = dhat::Profiler::new_heap();
    let curve = BezierCurve(Matrix2xX::from_iterator(4,
        [50.0, 200.0, 0.0, 50.0, 0.0, 33.0, 66.0, 100.0]
    ));
    black_box(curve.split(0.5));
    black_box(curve.castlejau_eval(0.5));
    black_box(curve);
    drop(profiler);
}