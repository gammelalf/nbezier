pub mod vector;
pub mod polynomial;
pub mod graham_scan;
pub mod bezier;

#[cfg(test)]
mod tests {
    use crate::bezier::{get_bernstein_polynomials, BezierCurve, pascal_triangle};
    use crate::graham_scan;
    use crate::polynomial::Polynomial;
    use crate::vector::Vector;

    #[test]
    fn bezier() {
        let line = BezierCurve(vec![Vector([0.0, 0.0]), Vector([1.0, 1.0])]);
        let (l, u) = line.split(0.5).unwrap();
        assert_eq!(l, BezierCurve(vec![Vector([0.0, 0.0]), Vector([0.5, 0.5])]));
        assert_eq!(u, BezierCurve(vec![Vector([0.5, 0.5]), Vector([1.0, 1.0])]));
    }

    #[test]
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
    }

    #[test]
    fn polynomials() {
        assert_eq!(Polynomial(vec![1.0, 2.0, 3.0]).derive(), Polynomial(vec![2.0, 6.0]));

        let cmp_float = |x: &f64, y: &f64| x.partial_cmp(y).expect("A wild NaN appeared");

        // Roots:
        for (n, m) in [
            (-3.0, 11.0),
            (-2.1, 3.9),
            (1.0, 2.0),
            (2.0, 3.0),
            (4.0, 5.0),
        ] {
            let p = &Polynomial(vec![-n, 1.0]) * &Polynomial(vec![-m, 1.0]);
            let mut roots = p.roots(); roots.sort_by(cmp_float);
            assert_eq!(roots, vec![n, m])
        }
    }

    #[test]
    fn pascal() {
        assert_eq!(pascal_triangle(0), vec![1]);
        assert_eq!(pascal_triangle(1), vec![1, 1]);
        assert_eq!(pascal_triangle(2), vec![1, 2, 1]);
        assert_eq!(pascal_triangle(3), vec![1, 3, 3, 1]);
        assert_eq!(pascal_triangle(4), vec![1, 4, 6, 4, 1]);
        assert_eq!(pascal_triangle(5), vec![1, 5, 10, 10, 5, 1]);
    }

    #[test]
    fn bernstein() {
        assert_eq!(get_bernstein_polynomials(0), vec![
            Polynomial(vec![1]),
        ]);
        assert_eq!(get_bernstein_polynomials(1), vec![
            Polynomial(vec![1, -1]),
            Polynomial(vec![0,  1]),
        ]);
        assert_eq!(get_bernstein_polynomials(2), vec![
            Polynomial(vec![1, -2,  1]),
            Polynomial(vec![0,  2, -2]),
            Polynomial(vec![0,  0,  1]),
        ]);
        assert_eq!(get_bernstein_polynomials(3), vec![
            Polynomial(vec![1, -3,  3, -1]),
            Polynomial(vec![0,  3, -6,  3]),
            Polynomial(vec![0,  0,  3, -3]),
            Polynomial(vec![0,  0,  0,  1]),
        ]);
    }
}
