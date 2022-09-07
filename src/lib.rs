#![warn(missing_docs)]
#![doc = include_str!("../README.md")]

pub mod bounding_box;
pub mod graham_scan;
pub mod nbezier;
pub mod npolynomial;
pub mod simple;

pub use crate::nbezier::BezierCurve;
pub use crate::simple::SimpleCurve;

#[cfg(test)]
mod tests {
    use crate::nbezier::{bernstein_polynomials, BezierCurve};
    use nalgebra::{
        Const, Matrix1, Matrix2, Matrix2x4, Matrix2x5, Matrix3, Matrix4, RowVector1, RowVector2,
        RowVector3, RowVector4, Vector2,
    };
    //use crate::graham_scan;
    use crate::npolynomial::Polynomial;

    #[test]
    fn bezier_split() {
        let line = BezierCurve(Matrix2::from_columns(&[
            Vector2::new(0.0, 0.0),
            Vector2::new(1.0, 1.0),
        ]));
        let [l, u] = line.split(0.5);
        assert_eq!(
            l.0,
            Matrix2::from_columns(&[Vector2::new(0.0, 0.0), Vector2::new(0.5, 0.5)])
        );
        assert_eq!(
            u.0,
            Matrix2::from_columns(&[Vector2::new(0.5, 0.5), Vector2::new(1.0, 1.0)])
        );
    }

    #[test]
    fn bezier_locate_point() {
        let curve = BezierCurve(Matrix2x4::from_columns(&[
            Vector2::new(50.0, 0.0),
            Vector2::new(200.0, 33.0),
            Vector2::new(0.0, 66.0),
            Vector2::new(50.0, 100.0),
        ]));
        for i in 1..10 {
            let actual_t = i as f64 / 10.0;
            let point = curve.castlejau_eval(actual_t);
            if let Some(found_t) = curve.locate_point(point) {
                assert!((found_t - actual_t).abs() < 1e-4);
            } else {
                panic!("{} not found!", actual_t);
            }
        }
    }

    #[test]
    fn bezier_derivative() {
        let curve = BezierCurve(Matrix2x5::from_columns(&[
            Vector2::new(50.0, 0.0),
            Vector2::new(200.0, 33.0),
            Vector2::new(0.0, 66.0),
            Vector2::new(200.0, 33.0),
            Vector2::new(50.0, 100.0),
        ]));
        assert_eq!(curve.derivative(), curve.polynomial().derive());
    }

    /*#[test]
    fn square() {
        // Counter clock wise
        let points: Vec<Vector<i8, 2>> = vec![Vector([0,0]), Vector([1,0]), Vector([1,1]), Vector([0,1])];
        let hull = graham_scan::convex_hull(points.clone());
        assert_eq!(points, hull);

        // Clock wise
        let points: Vec<Vector<f32, 2>> = vec![Vector([0.0, 0.0]), Vector([0.0, 1.0]), Vector([1.0, 1.0]), Vector([1.0, 0.0])];
        let hull = graham_scan::convex_hull(points.clone());
        assert_ne!(points, hull);
    }

    #[test]
    fn six_points() {
        let points: Vec<Vector<i32, 2>> = vec![
            Vector([ 6,  2]),
            Vector([ 4,  7]),
            Vector([10,  5]),
            Vector([12,  2]),
            Vector([10, 10]),
            Vector([15,  7]),
        ];
        let hull = graham_scan::convex_hull(points);
        assert_eq!(hull, vec![
            Vector([ 6,  2]),
            Vector([12,  2]),
            Vector([15,  7]),
            Vector([10, 10]),
            Vector([ 4,  7]),
        ]);
    }*/

    #[test]
    fn polynomials() {
        assert_eq!(
            Polynomial(RowVector3::new(1.0, 2.0, 3.0)).derive(),
            Polynomial(RowVector2::new(2.0, 6.0))
        );

        let cmp_float = |x: &f64, y: &f64| x.partial_cmp(y).expect("A wild NaN appeared");

        // Roots:
        for (n, m) in [
            (-3.0, 11.0),
            (-2.1, 3.9),
            (1.0, 2.0),
            (2.0, 3.0),
            (4.0, 5.0),
        ] {
            let p = Polynomial(RowVector2::new(-n, 1.0)).mul(&RowVector2::new(-m, 1.0));
            let mut roots = p.roots();
            roots.sort_by(cmp_float);
            assert_eq!(roots, vec![n, m])
        }
    }

    #[test]
    fn bernstein() {
        assert_eq!(
            bernstein_polynomials::<f32, _>(Const::<1>),
            Matrix1::from_rows(&[RowVector1::new(1.0),])
        );
        assert_eq!(
            bernstein_polynomials::<f32, _>(Const::<2>),
            Matrix2::from_rows(&[RowVector2::new(1.0, -1.0), RowVector2::new(0.0, 1.0),])
        );
        assert_eq!(
            bernstein_polynomials::<f32, _>(Const::<3>),
            Matrix3::from_rows(&[
                RowVector3::new(1.0, -2.0, 1.0),
                RowVector3::new(0.0, 2.0, -2.0),
                RowVector3::new(0.0, 0.0, 1.0),
            ])
        );
        assert_eq!(
            bernstein_polynomials::<f32, _>(Const::<4>),
            Matrix4::from_rows(&[
                RowVector4::new(1.0, -3.0, 3.0, -1.0),
                RowVector4::new(0.0, 3.0, -6.0, 3.0),
                RowVector4::new(0.0, 0.0, 3.0, -3.0),
                RowVector4::new(0.0, 0.0, 0.0, 1.0),
            ])
        );
    }
}
