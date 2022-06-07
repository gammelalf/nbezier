use once_cell::sync::Lazy;
use smallvec::smallvec;
use nalgebra::Vector2;
use gammalg::bezier::BezierCurve;

pub static CURVES: Lazy<Vec<BezierCurve<f64>>> = Lazy::new(|| vec![
    BezierCurve(smallvec![
        Vector2::new(0.0, 0.0),
        Vector2::new(10.0, 0.0),
        Vector2::new(0.0, 10.0),
        Vector2::new(10.0, 10.0),
    ]),
    BezierCurve(smallvec![
        Vector2::new(50.0, 0.0),
        Vector2::new(200.0, 33.0),
        Vector2::new(0.0, 66.0),
        Vector2::new(50.0, 100.0),
    ]),
    BezierCurve(smallvec![
        Vector2::new(50.0, 0.0),
        Vector2::new(155.0, 23.1),
        Vector2::new(88.5, 46.2),
        Vector2::new(56.3, 69.643)
    ]),
    BezierCurve(smallvec![
        Vector2::new(56.3, 69.643),
        Vector2::new(42.5, 79.69),
        Vector2::new(35.0, 89.8),
        Vector2::new(50.0, 100.0)
    ]),
]);
