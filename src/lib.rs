pub mod bezier;
pub mod bounding_box;
pub mod graham_scan;
pub mod npolynomial;

pub use bezier::BezierCurve;

#[cfg(test)]
mod tests {
    use crate::bezier::{bernstein_polynomials, pascal_triangle, BezierCurve};
    use nalgebra::{RowDVector, RowVector2, RowVector3, Vector2};
    use smallvec::smallvec;
    //use crate::graham_scan;
    use crate::npolynomial::Polynomial;

    #[test]
    fn bezier_split() {
        let line = BezierCurve(smallvec![Vector2::new(0.0, 0.0), Vector2::new(1.0, 1.0)]);
        let (l, u) = line.split(0.5).unwrap();
        assert_eq!(
            l,
            BezierCurve(smallvec![Vector2::new(0.0, 0.0), Vector2::new(0.5, 0.5)])
        );
        assert_eq!(
            u,
            BezierCurve(smallvec![Vector2::new(0.5, 0.5), Vector2::new(1.0, 1.0)])
        );
    }

    #[test]
    fn bezier_locate_point() {
        let curve = BezierCurve(smallvec![
            Vector2::new(50.0, 0.0),
            Vector2::new(200.0, 33.0),
            Vector2::new(0.0, 66.0),
            Vector2::new(50.0, 100.0),
        ]);
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
        let curve = BezierCurve(smallvec![
            Vector2::new(50.0, 0.0),
            Vector2::new(200.0, 33.0),
            Vector2::new(0.0, 66.0),
            Vector2::new(50.0, 100.0),
        ]);
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
    fn pascal() {
        assert_eq!(pascal_triangle::<i32>(0), vec![1]);
        assert_eq!(pascal_triangle::<i32>(1), vec![1, 1]);
        assert_eq!(pascal_triangle::<i32>(2), vec![1, 2, 1]);
        assert_eq!(pascal_triangle::<i32>(3), vec![1, 3, 3, 1]);
        assert_eq!(pascal_triangle::<i32>(4), vec![1, 4, 6, 4, 1]);
        assert_eq!(pascal_triangle::<i32>(5), vec![1, 5, 10, 10, 5, 1]);
    }

    #[test]
    fn bernstein() {
        assert_eq!(
            bernstein_polynomials::<f32>(0),
            vec![Polynomial(RowDVector::from_row_slice(&[1.0])),]
        );
        assert_eq!(
            bernstein_polynomials::<f32>(1),
            vec![
                Polynomial(RowDVector::from_row_slice(&[1.0, -1.0])),
                Polynomial(RowDVector::from_row_slice(&[0.0, 1.0])),
            ]
        );
        assert_eq!(
            bernstein_polynomials::<f32>(2),
            vec![
                Polynomial(RowDVector::from_row_slice(&[1.0, -2.0, 1.0])),
                Polynomial(RowDVector::from_row_slice(&[0.0, 2.0, -2.0])),
                Polynomial(RowDVector::from_row_slice(&[0.0, 0.0, 1.0])),
            ]
        );
        assert_eq!(
            bernstein_polynomials::<f32>(3),
            vec![
                Polynomial(RowDVector::from_row_slice(&[1.0, -3.0, 3.0, -1.0])),
                Polynomial(RowDVector::from_row_slice(&[0.0, 3.0, -6.0, 3.0])),
                Polynomial(RowDVector::from_row_slice(&[0.0, 0.0, 3.0, -3.0])),
                Polynomial(RowDVector::from_row_slice(&[0.0, 0.0, 0.0, 1.0])),
            ]
        );
    }
}
