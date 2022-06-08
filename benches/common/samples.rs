use once_cell::sync::Lazy;
use smallvec::smallvec;
use nalgebra::Vector2;
use gammalg::bezier::BezierCurve;

/// Points generated randomly
/// ```python
/// from random import random
/// for i in range(10):
///     print(f"Vector2::new({(random()-0.5)*i}, {(random()-0.5)*i})")
/// ```
pub static POINTS: [Vector2<f64>; 10] = [
    Vector2::new( 0.0,      0.0    ),
    Vector2::new(-0.29734,  0.44984),
    Vector2::new(-0.52560,  0.42885),
    Vector2::new( 1.42777, -0.02652),
    Vector2::new( 1.98032, -0.67824),
    Vector2::new( 0.44863, -0.91328),
    Vector2::new(-2.51139, -0.79100),
    Vector2::new(-3.10479, -0.59318),
    Vector2::new(-1.16022, -2.95591),
    Vector2::new(-1.07946,  0.78888),
];

pub static CURVES: Lazy<Curves> = Lazy::new(Curves::new);
#[allow(non_snake_case)]
pub struct Curves {
    pub LINEAR: Vec<BezierCurve<f64>>,
    pub QUADRATIC: Vec<BezierCurve<f64>>,
    pub CUBIC: Vec<BezierCurve<f64>>,
    pub HIGHER: Vec<BezierCurve<f64>>,
}
impl Curves {
    pub fn new() -> Curves {
        Curves {
            LINEAR:    vec![
                BezierCurve([0, 1].into_iter().map(|i| POINTS[i]).collect()),
                BezierCurve([1, 2].into_iter().map(|i| POINTS[i]).collect()),
                BezierCurve([2, 3].into_iter().map(|i| POINTS[i]).collect()),
                BezierCurve([3, 4].into_iter().map(|i| POINTS[i]).collect()),
                BezierCurve([4, 5].into_iter().map(|i| POINTS[i]).collect()),
                BezierCurve([5, 6].into_iter().map(|i| POINTS[i]).collect()),
                BezierCurve([6, 7].into_iter().map(|i| POINTS[i]).collect()),
                BezierCurve([7, 8].into_iter().map(|i| POINTS[i]).collect()),
                BezierCurve([8, 9].into_iter().map(|i| POINTS[i]).collect()),
                BezierCurve([9, 0].into_iter().map(|i| POINTS[i]).collect()),
            ],
            QUADRATIC: vec![
                BezierCurve([0, 1, 2].into_iter().map(|i| POINTS[i]).collect()),
                BezierCurve([2, 3, 4].into_iter().map(|i| POINTS[i]).collect()),
                BezierCurve([4, 6, 8].into_iter().map(|i| POINTS[i]).collect()),
                BezierCurve([1, 3, 5].into_iter().map(|i| POINTS[i]).collect()),
                BezierCurve([5, 7, 9].into_iter().map(|i| POINTS[i]).collect()),
                BezierCurve([0, 2, 1].into_iter().map(|i| POINTS[i]).collect()),
                BezierCurve([2, 4, 3].into_iter().map(|i| POINTS[i]).collect()),
                BezierCurve([4, 8, 6].into_iter().map(|i| POINTS[i]).collect()),
                BezierCurve([1, 5, 3].into_iter().map(|i| POINTS[i]).collect()),
                BezierCurve([5, 9, 7].into_iter().map(|i| POINTS[i]).collect()),
            ],
            CUBIC:     vec![
                BezierCurve([0, 1, 2, 3].into_iter().map(|i| POINTS[i]).collect()),
                BezierCurve([2, 3, 4, 5].into_iter().map(|i| POINTS[i]).collect()),
                BezierCurve([4, 5, 6, 7].into_iter().map(|i| POINTS[i]).collect()),
                BezierCurve([6, 7, 8, 9].into_iter().map(|i| POINTS[i]).collect()),
                BezierCurve([8, 9, 0, 1].into_iter().map(|i| POINTS[i]).collect()),
                BezierCurve([0, 4, 2, 6].into_iter().map(|i| POINTS[i]).collect()),
                BezierCurve([1, 5, 3, 7].into_iter().map(|i| POINTS[i]).collect()),
                BezierCurve([2, 6, 4, 8].into_iter().map(|i| POINTS[i]).collect()),
                BezierCurve([3, 7, 5, 9].into_iter().map(|i| POINTS[i]).collect()),
                BezierCurve([4, 8, 6, 0].into_iter().map(|i| POINTS[i]).collect()),
            ],
            HIGHER:    vec![
                BezierCurve([0, 1, 2, 3, 4].into_iter().map(|i| POINTS[i]).collect()),
                BezierCurve([2, 3, 4, 5, 6].into_iter().map(|i| POINTS[i]).collect()),
                BezierCurve([4, 5, 6, 7, 8].into_iter().map(|i| POINTS[i]).collect()),
                BezierCurve([6, 7, 8, 9, 0].into_iter().map(|i| POINTS[i]).collect()),
                BezierCurve([8, 9, 0, 1, 2].into_iter().map(|i| POINTS[i]).collect()),
                BezierCurve([0, 4, 2, 6, 8].into_iter().map(|i| POINTS[i]).collect()),
                BezierCurve([1, 5, 3, 7, 9].into_iter().map(|i| POINTS[i]).collect()),
                BezierCurve([2, 6, 4, 8, 0].into_iter().map(|i| POINTS[i]).collect()),
                BezierCurve([3, 7, 5, 9, 1].into_iter().map(|i| POINTS[i]).collect()),
                BezierCurve([4, 8, 6, 0, 2].into_iter().map(|i| POINTS[i]).collect()),
            ],
        }
    }

    pub fn iter(&self) -> impl Iterator<Item=&BezierCurve<f64>> {
        self.LINEAR.iter()
            .chain(self.QUADRATIC.iter())
            .chain(self.CUBIC.iter())
            .chain(self.HIGHER.iter())
    }
}
