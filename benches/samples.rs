use gammalg::bezier::BezierCurve;
use gammalg::vector::Vector;
use once_cell::sync::Lazy;

pub static CURVES: Lazy<Vec<BezierCurve<f64>>> = Lazy::new(|| vec![
    BezierCurve(vec![
        Vector([0.0, 0.0]),
        Vector([10.0, 0.0]),
        Vector([0.0, 10.0]),
        Vector([10.0, 10.0]),
    ]),
    BezierCurve(vec![
        Vector([50.0, 0.0]),
        Vector([200.0, 33.0]),
        Vector([0.0, 66.0]),
        Vector([50.0, 100.0]),
    ]),
    BezierCurve(vec![
        Vector([50.0, 0.0]),
        Vector([155.0, 23.1]),
        Vector([88.5, 46.2]),
        Vector([56.3, 69.643])
    ]),
    BezierCurve(vec![
        Vector([56.3, 69.643]),
        Vector([42.5, 79.69]),
        Vector([35.0, 89.8]),
        Vector([50.0, 100.0])
    ]),
]);
