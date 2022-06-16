use crate::bounding_box::BoundingBox;
use crate::graham_scan::convex_hull;
use crate::npolynomial::{Polynomial, Polynomial1xX, Polynomial2xX};
use nalgebra::{
    ComplexField, Dynamic, Field, Matrix2xX, OMatrix, RealField, RowDVector, RowVector2,
    RowVector3, RowVector4, Scalar, Vector2,
};
use num::{Num, One};
use smallvec::{smallvec, SmallVec};
use std::ops::{Add, Deref, DerefMut};

#[derive(Clone, Debug, PartialEq)]
pub struct BezierCurve<T: Scalar>(pub CurveInternal<T>);
type CurveInternal<T> = SmallVec<[Vector2<T>; 4]>;

impl<T: Scalar> Deref for BezierCurve<T> {
    type Target = CurveInternal<T>;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}
impl<T: Scalar> DerefMut for BezierCurve<T> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl<T: Scalar> BezierCurve<T> {
    /// Returns a curve's degree which is one lower then its number of control points
    pub fn degree(&self) -> usize {
        self.len() - 1
    }
}

impl<T: ComplexField + Scalar + std::fmt::Display> BezierCurve<T> {
    /// Returns a new bezier curve of the same shape whose degree is one step higher
    pub fn raise_degree(&self) -> BezierCurve<T> {
        // Transform points
        let old = Matrix2xX::from_columns(&self[..]).transpose();
        let new = BezierCurve::elevation_matrix(self.len()) * &old;

        // Create new SmallVec from matrix
        let mut points = SmallVec::with_capacity(self.len() - 1);
        for row in new.row_iter() {
            points.push(row.transpose());
        }
        BezierCurve(points)
    }

    /// Returns a new bezier curve of an approximated shape whose degree is one step lower
    pub fn lower_degree(&self) -> BezierCurve<T> {
        // Compute (M^T M)^{-1} M^T
        let matrix = BezierCurve::elevation_matrix(self.len() - 1);
        let transposed = matrix.transpose();
        let square = matrix.tr_mul(&matrix);
        let inverse = square.try_inverse().unwrap();
        let matrix = inverse * transposed;

        // Transform points
        let old = Matrix2xX::from_columns(&self[..]).transpose();
        let new = &matrix * &old;

        // Create new SmallVec from matrix
        let mut points = SmallVec::with_capacity(self.len() - 1);
        for row in new.row_iter() {
            points.push(row.transpose());
        }
        BezierCurve(points)
    }

    /// Constructs the matrix to raise a curve's degree from n to n+1
    fn elevation_matrix(n: usize) -> OMatrix<T, Dynamic, Dynamic> {
        let n_plus_one = convert::usize_to_generic::<T>(n + 1);
        let mut matrix = OMatrix::zeros_generic(Dynamic::new(n + 1), Dynamic::new(n));
        matrix[(0, 0)] = T::one();
        matrix[(n, n - 1)] = T::one();
        let mut i_as_k = T::zero();
        for i in 1..n {
            i_as_k += T::one();
            matrix[(i, i - 1)] = i_as_k.clone() / n_plus_one.clone();
            matrix[(i, i)] = T::one() - matrix[(i, i - 1)].clone();
        }
        matrix
    }
}

/* Hitboxes and intersections */
impl<T: RealField + Scalar> BezierCurve<T> {
    /// Constructs an axis aligned bounding box containing all controll points.
    ///
    /// This box will also contain the whole curve, but can highly overestimate it.
    /// It can be used as a fast way to estimate intersections.
    /// For more precise checks consider: `minimal_bounding_box` or `convex_hull`
    pub fn bounding_box(&self) -> BoundingBox<T> {
        BoundingBox::from_slice(&self)
    }

    /// Computes the smallest possible axis aligend bounding box containing the curve.
    ///
    /// As this relies on computing the derivative's roots, only cubic curves and lower are
    /// implemented so far.
    /// **Higher curves than cubics will panic!**
    pub fn minimal_bounding_box(&self) -> BoundingBox<T> {
        assert!(self.degree() < 4);
        let mut points: SmallVec<[Vector2<T>; 6]> = SmallVec::new();
        points.push(self[0].clone());
        points.push(self[self.len() - 1].clone());
        let d = self.derivative();
        let roots_x = Polynomial(d.0.row(0)).roots();
        let roots_y = Polynomial(d.0.row(1)).roots();
        roots_x
            .into_iter()
            .chain(roots_y.into_iter())
            .filter(|t| &T::zero() <= t && t <= &T::one())
            .for_each(|t| {
                points.push(self.castlejau_eval(t));
            });
        BoundingBox::from_iter(points.into_iter())
    }

    /// Computes the control points' convex hull using [graham scan](https://en.wikipedia.org/wiki/Graham_scan).
    ///
    /// *Actual code for intersecting convex polygons is not implemented yet!*
    pub fn convex_hull(&self) -> Vec<Vector2<T>> {
        convex_hull(self.0.clone().into_vec())
    }

    /// Trys to locate a point on the curve and returns its `t` value if located.
    ///
    /// This function was mostly ported from [this](https://github.com/dhermes/bezier) python
    /// package.
    ///
    /// It repeatedly subdivides the curve discarding subcurves whose `bounding_box` doesn't
    /// contain the point.
    /// After a set number of subdivisions, it collects all remaining subcurves into an
    /// approximation for `t`.
    /// Finally it refines this approximation by running newton's method a set number of iterators.
    ///
    /// Currently it runs 20 subdivisions and 2 iterations of newton.
    pub fn locate_point(&self, p: Vector2<T>) -> Option<T> {
        let zero = T::zero();
        let one = T::one();
        let two = one.clone() + one.clone();
        let halve = one.clone() / two.clone();

        #[allow(non_snake_case)]
        let NUM_SUBDIVISIONS: usize = 20;
        #[allow(non_snake_case)]
        let MAX_DEVIATION: T = halve.powi(19);
        #[allow(non_snake_case)]
        let NUM_NEWTON: usize = 2;

        // Check basic bounding box before doing any allocations
        let bb = self.bounding_box();
        if !bb.contains(p.clone()) {
            return None;
        }

        // Subdivide curve and check bounding boxes
        let mut divisions = (&mut vec![SubCurve::from(self.clone())], &mut Vec::new());
        for _ in 0..NUM_SUBDIVISIONS {
            if divisions.0.len() == 0 {
                return None;
            }
            for curve in divisions.0.iter() {
                let bb = curve.bounding_box();
                if bb.contains(p.clone()) {
                    let (lower, upper) = curve.split();
                    divisions.1.push(lower);
                    divisions.1.push(upper);
                }
            }
            divisions.0.clear();
            divisions = (divisions.1, divisions.0);
        }
        let divisions = divisions.0;
        if divisions.len() == 0 {
            return None;
        }

        // Combine all subdivisions' start and end into a single approximation
        let mut divisions_len = zero.clone(); // Replacement for divisions.len() which would be usize
        let mut approximation = zero.clone();
        for curve in divisions.iter() {
            approximation = approximation + curve.from.clone() + curve.to.clone();
            divisions_len = divisions_len + one.clone();
        }
        approximation = approximation.clone() / (two.clone() * divisions_len.clone());

        // Check deviation to see if subdivisions are too far apart
        let mut deviation = zero.clone();
        for curve in divisions.iter() {
            deviation = deviation.clone()
                + (curve.from.clone() - approximation.clone()).powi(2)
                + (curve.to.clone() - approximation.clone()).powi(2);
        }
        deviation = (deviation / (two * divisions_len.clone())).sqrt();
        if deviation > MAX_DEVIATION {
            return None;
        }

        // Refine solution using newton's method
        for _ in 0..NUM_NEWTON {
            let function = self.castlejau_eval(approximation.clone());
            let derivative = self.tangent(approximation.clone());
            let f_minus_p: Vector2<T> = function.clone() - p.clone();
            let d_times_d: T = derivative.dot(&derivative);
            let d_times_f_minus_p: T = derivative.dot(&f_minus_p);
            approximation = approximation - d_times_f_minus_p / d_times_d;
            //approximation = approximation - (derivative * (function - p)) / (derivative * derivative);
        }

        Some(approximation)
    }

    /// WIP
    pub fn get_intersections(&self, other: &Self) -> Vec<Vector2<T>> {
        #[allow(non_snake_case)]
        let NUM_SUBDIVISIONS: usize = 20;

        // Subdivide curves and check bounding boxes
        let mut divisions = (
            &mut vec![(SubCurve::from(self.clone()), SubCurve::from(other.clone()))],
            &mut Vec::new(),
        );
        for _ in 0..NUM_SUBDIVISIONS {
            if divisions.0.len() == 0 {
                return Vec::new();
            }
            for (s, o) in divisions.0.iter() {
                if s.bounding_box().intersects(&o.bounding_box()) {
                    let (s_lower, s_upper) = s.split();
                    let (o_lower, o_upper) = o.split();
                    divisions.1.push((s_lower.clone(), o_lower.clone()));
                    divisions.1.push((s_lower, o_upper.clone()));
                    divisions.1.push((s_upper.clone(), o_lower));
                    divisions.1.push((s_upper, o_upper));
                }
            }
            divisions.0.clear();
            divisions = (divisions.1, divisions.0);
        }
        let divisions = divisions.0;
        if divisions.len() == 0 {
            return Vec::new();
        }

        return divisions.iter().map(|(s, _)| s.middle_point()).collect();
    }
}

/* Stuff using de castlejau 's algorithm */
impl<T: Field + Scalar> BezierCurve<T> {
    /// Splits a curve into two parts
    ///
    /// The first part is the same shape as the original curve between 0 and t and the second
    /// part as the curve between t and 1.
    /// This method assumes `t` to between 0 and 1 but doesn't check it.
    pub fn split(&self, t: T) -> (BezierCurve<T>, BezierCurve<T>) {
        let inv_t = T::one() - t.clone();
        match &self[..] {
            [] | [_] => (self.clone(), self.clone()),
            [a2, b2] => {
                let a1 = a2 * inv_t + b2 * t;
                (
                    BezierCurve(smallvec![a2.clone(), a1.clone()]),
                    BezierCurve(smallvec![a1, b2.clone()]),
                )
            }
            [a3, b3, c3] => {
                let a2 = a3 * inv_t.clone() + b3 * t.clone();
                let b2 = b3 * inv_t.clone() + c3 * t.clone();
                let a1 = &a2 * inv_t + &b2 * t;
                (
                    BezierCurve(smallvec![a3.clone(), a2, a1.clone()]),
                    BezierCurve(smallvec![a1, b2, c3.clone()]),
                )
            }
            [a4, b4, c4, d4] => {
                let a3 = a4 * inv_t.clone() + b4 * t.clone();
                let b3 = b4 * inv_t.clone() + c4 * t.clone();
                let c3 = c4 * inv_t.clone() + d4 * t.clone();
                let a2 = &a3 * inv_t.clone() + &b3 * t.clone();
                let b2 = &b3 * inv_t.clone() + &c3 * t.clone();
                let a1 = &a2 * inv_t + &b2 * t;
                (
                    BezierCurve(smallvec![a4.clone(), a3, a2, a1.clone()]),
                    BezierCurve(smallvec![a1, b2, c3, d4.clone()]),
                )
            }
            _ => {
                let len = self.len();

                let mut lower = SmallVec::with_capacity(len);
                let mut upper = SmallVec::with_capacity(len);
                lower.push(self[0].clone());
                upper.push(self[len - 1].clone());

                let mut points = (&mut self.0.clone(), &mut self.0.clone());
                for i in 1..len {
                    BezierCurve::castlejau_step(points.0, points.1, t.clone());
                    lower.push(points.1[0].clone());
                    upper.push(points.1[len - i - 1].clone());
                    points = (points.1, points.0);
                }
                upper.reverse(); // I find it more intuitive if the t goes through the two parts in the same direction
                (BezierCurve(lower), BezierCurve(upper))
            }
        }
    }

    /// Get the point on the curve at position `t`.
    ///
    /// This method uses de castlejau's algorithm. An alternative way would be to evaluate the
    /// curve's polynomial (See `BezierCurve::polynomial`).
    pub fn castlejau_eval(&self, t: T) -> Vector2<T> {
        let inv_t = T::one() - t.clone();
        match &self[..] {
            [] => panic!(),
            [a1] => a1.clone(),
            [a2, b2] => a2 * inv_t + b2 * t,
            [a3, b3, c3] => {
                let a2 = a3 * inv_t.clone() + b3 * t.clone();
                let b2 = b3 * inv_t.clone() + c3 * t.clone();
                a2 * inv_t + b2 * t
            }
            [a4, b4, c4, d4] => {
                let a3 = a4 * inv_t.clone() + b4 * t.clone();
                let b3 = b4 * inv_t.clone() + c4 * t.clone();
                let c3 = c4 * inv_t.clone() + d4 * t.clone();
                let a2 = &a3 * inv_t.clone() + &b3 * t.clone();
                let b2 = &b3 * inv_t.clone() + &c3 * t.clone();
                a2 * inv_t + b2 * t
            }
            _ => {
                let mut old_points = self.0.clone();
                let mut new_points = self.0.clone();
                let mut points = (&mut old_points, &mut new_points);
                while points.1.len() > 1 {
                    BezierCurve::castlejau_step(points.0, points.1, t.clone());
                    points = (points.1, points.0);
                }
                return points.1[0].clone();
            }
        }
    }

    /// Performs a single step of de castlejau's algorithm
    ///
    /// i.e. combines `n` points into `n - 1` points by computing `(1 - t) * A + t * B` on
    /// consecutive points `A` and `B`
    fn castlejau_step(input: &CurveInternal<T>, output: &mut CurveInternal<T>, t: T) {
        output.clear();
        let len = input.len();
        let t_inv = T::one() - t.clone();
        for (p, q) in (&input[0..len - 1]).iter().zip((&input[1..len]).iter()) {
            output.push(p * t_inv.clone() + q * t.clone());
        }
    }
}

/* Stuff using polynomials */
impl<T: Field + Scalar> BezierCurve<T> {
    /// Computes the curve's polynomial
    ///
    /// This polynomial evaluated between 0 and 1 yields the same points as its corrisponding bezier curve.
    pub fn polynomial(&self) -> Polynomial2xX<T> {
        let zero = T::zero();
        let one = T::one();
        match &self[..] {
            [a] => Polynomial(convert::static_to_dynamic(a)),
            [a, b] => {
                let p_a = a * RowVector2::new(one.clone(), zero.clone() - one.clone());
                let p_b = b * RowVector2::new(zero, one);
                let p = p_a + p_b;
                Polynomial(convert::static_to_dynamic(&p))
            }
            [a, b, c] => {
                let two = one.clone() + one.clone();
                let p_a = a * RowVector3::new(one.clone(), zero.clone() - two.clone(), one.clone());
                let p_b = b * RowVector3::new(zero.clone(), two.clone(), zero.clone() - two);
                let p_c = c * RowVector3::new(zero.clone(), zero, one);
                let p = p_a + p_b + p_c;
                Polynomial(convert::static_to_dynamic(&p))
            }
            [a, b, c, d] => {
                let three = one.clone() + one.clone() + one.clone();
                let six = three.clone() + three.clone();
                let p_a = a * RowVector4::new(
                    one.clone(),
                    zero.clone() - three.clone(),
                    three.clone(),
                    zero.clone() - one.clone(),
                );
                let p_b = b * RowVector4::new(
                    zero.clone(),
                    three.clone(),
                    zero.clone() - six,
                    three.clone(),
                );
                let p_c = c * RowVector4::new(
                    zero.clone(),
                    zero.clone(),
                    three.clone(),
                    zero.clone() - three,
                );
                let p_d = d * RowVector4::new(zero.clone(), zero.clone(), zero, one);
                let p = p_a + p_b + p_c + p_d;
                Polynomial(convert::static_to_dynamic(&p))
            }
            _ => {
                let mut ps = bernstein_polynomials::<T>(self.degree())
                    .into_iter()
                    .zip(self.iter())
                    .map(|(p, a)| a * p.0);
                if let Some(mut p) = ps.next() {
                    for q in ps {
                        p += q;
                    }
                    return Polynomial(p);
                } else {
                    unreachable!();
                };
            }
        }
    }

    /// Computes the curve's polynomial's derivative
    ///
    /// This method is a faster alternative to calling `Polynomial::derive` on the result of
    /// `BezierCurve::polynomial`.
    pub fn derivative(&self) -> Polynomial2xX<T> {
        let zero = T::zero();
        let one = T::one();
        match &self[..] {
            [_] => Polynomial(Matrix2xX::zeros(0)),
            [a, b] => {
                let p = b - a;
                Polynomial(convert::static_to_dynamic(&p))
            }
            [a, b, c] => {
                let two = one.clone() + one.clone();
                let p_a = (b - a) * RowVector2::new(one.clone(), zero.clone() - one.clone());
                let p_b = (c - b) * RowVector2::new(zero, one);
                let p = (p_a + p_b) * two;
                Polynomial(convert::static_to_dynamic(&p))
            }
            [a, b, c, d] => {
                let two = one.clone() + one.clone();
                let three = two.clone() + one.clone();
                let p_a =
                    (b - a) * RowVector3::new(one.clone(), zero.clone() - two.clone(), one.clone());
                let p_b = (c - b) * RowVector3::new(zero.clone(), two.clone(), zero.clone() - two);
                let p_c = (d - c) * RowVector3::new(zero.clone(), zero, one);
                let p = (p_a + p_b + p_c) * three;
                Polynomial(convert::static_to_dynamic(&p))
            }
            _ => {
                let mut ps = bernstein_polynomials::<T>(self.degree() - 1)
                    .into_iter()
                    .enumerate()
                    .map(|(i, p)| {
                        let a = &self[i + 1] - &self[i];
                        a * p.0
                    });
                if let Some(mut p) = ps.next() {
                    for q in ps {
                        p += q;
                    }
                    return Polynomial(p);
                } else {
                    unreachable!();
                }
            }
        }
    }

    /// Computes the curve's tangent vector at `t`
    ///
    /// If you need a lot of them at once, it is far more efficent to call
    /// `BezierCurve::derivative` once yourself and evaluate it multiple times,
    /// as this function would recompute the derivative for every vector.
    ///
    /// *The resulting vector is not normalized!*
    pub fn tangent(&self, t: T) -> Vector2<T> {
        self.derivative().evaluate(t.clone())
    }

    /// Computes the curve's normal vector at `t`
    ///
    /// Similarly to `BezierCurve::tangent` calling this multiple times to get a lot of vectors
    /// would be inefficent. Compute their tangent vectors (see `BezierCurve::tangent`) and rotate
    /// them yourself instead.
    ///
    /// *The resulting vector is not normalized!*
    pub fn normal(&self, t: T) -> Vector2<T> {
        let [[x, y]] = self.tangent(t).data.0;
        Vector2::new(T::zero() - y, x)
    }
}
mod convert {
    use nalgebra::{ArrayStorage, Const, Matrix, Matrix2xX, U2};
    use num::Num;

    /// Convert a static polynomial matrix into a dynamic one
    pub(crate) fn static_to_dynamic<T: Clone, const C: usize>(
        matrix: &Matrix<T, U2, Const<C>, ArrayStorage<T, 2, C>>,
    ) -> Matrix2xX<T> {
        let mut data = Vec::with_capacity(2 * C);
        data.extend_from_slice(matrix.data.as_slice());
        let data = nalgebra::VecStorage::new(nalgebra::Const::<2>, nalgebra::Dynamic::new(C), data);
        Matrix2xX::from_data(data)
    }

    /// Helper function used when a formula treats a curve's degree as a scalar
    pub(crate) fn usize_to_generic<T: Num>(n: usize) -> T {
        let mut k = T::zero();
        for _ in 0..n {
            k = k + T::one();
        }
        k
    }
}

/// Computes the bernstein polynomial basis for a given degree
pub fn bernstein_polynomials<N>(degree: usize) -> Vec<Polynomial1xX<N>>
where
    N: Field + Scalar,
{
    let one = N::one();
    let zero = N::zero();

    // Prepare the powers of x and (1-x)
    let mut powers = (
        Vec::with_capacity(degree + 1),
        Vec::with_capacity(degree + 1),
    );

    powers
        .0
        .push(Polynomial(RowDVector::from_row_slice(&[one.clone()]))); // TODO don't alloc p(x) = 1
    powers
        .1
        .push(Polynomial(RowDVector::from_row_slice(&[one.clone()])));
    powers.0.push(Polynomial(RowDVector::from_row_slice(&[
        zero.clone(),
        one.clone(),
    ])));
    powers.1.push(Polynomial(RowDVector::from_row_slice(&[
        one.clone(),
        zero - one,
    ])));

    for i in 1..degree {
        powers.0.push(powers.0[i].mul(&powers.0[1].0));
        powers.1.push(powers.1[i].mul(&powers.1[1].0));
    }

    // Combine powers into Bernstein polynomials
    let mut base = Vec::with_capacity(degree + 1);
    let pascal = pascal_triangle::<N>(degree);
    for (i, coeff) in pascal.into_iter().enumerate() {
        let mut b = powers.0[i].mul(&powers.1[degree - i].0);
        b.0 *= coeff; // Use MutAssign to multiply inplace and avoid reallocation
        base.push(b);
    }

    return base;
}

/// Computes a given layer of pascal's triangle
///
/// This function is used in `bernstein_polynomials` to compute all binomial coefficient for a
/// single `n`.
pub fn pascal_triangle<N>(layer: usize) -> Vec<N>
where
    N: Add<Output = N> + One + Clone,
{
    let one = N::one();
    let mut old_layer = Vec::with_capacity(layer + 1);
    let mut new_layer = Vec::with_capacity(layer + 1);
    new_layer.push(N::one());

    while new_layer.len() < new_layer.capacity() {
        old_layer.clone_from(&new_layer);

        new_layer.push(one.clone());
        let get = |i| old_layer.get(i as usize).map(|n| n).unwrap_or(&one).clone();
        for i in 1..new_layer.len() - 1 {
            new_layer[i] = get(i - 1) + get(i);
        }
    }

    new_layer
}

/* Helper struct repeatedly subdividing a curve */
#[derive(Clone)]
struct SubCurve<T: RealField> {
    from: T,
    to: T,
    curve: BezierCurve<T>,
}
impl<T: RealField> SubCurve<T> {
    fn split(&self) -> (SubCurve<T>, SubCurve<T>) {
        let two = T::one() + T::one();
        let middle = (self.from.clone() + self.to.clone()) / two.clone();

        let (lower, upper) = self.curve.split(T::one() / two);
        (
            SubCurve {
                from: self.from.clone(),
                to: middle.clone(),
                curve: lower,
            },
            SubCurve {
                from: middle,
                to: self.to.clone(),
                curve: upper,
            },
        )
    }

    fn middle_point(&self) -> Vector2<T> {
        let two = T::one() + T::one();
        let middle = (self.from.clone() + self.to.clone()) / two;

        self.curve.castlejau_eval(middle)
    }
}
impl<T: RealField> From<BezierCurve<T>> for SubCurve<T> {
    fn from(curve: BezierCurve<T>) -> Self {
        SubCurve {
            from: T::zero(),
            to: T::one(),
            curve,
        }
    }
}
impl<T: RealField> Deref for SubCurve<T> {
    type Target = BezierCurve<T>;
    fn deref(&self) -> &Self::Target {
        &self.curve
    }
}
impl<T: RealField> DerefMut for SubCurve<T> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.curve
    }
}
