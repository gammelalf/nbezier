use crate::bounding_box::BoundingBox;
use crate::graham_scan::convex_hull;
use crate::npolynomial::{Polynomial, Polynomial1xX, Polynomial2xX};
use nalgebra::{
    ComplexField, Dynamic, Field, Matrix, Matrix2xX, RealField, RowDVector, RowVector2, RowVector3,
    RowVector4, Scalar, Vector2,
};
use num::{Num, One};
use smallvec::{smallvec, SmallVec};
use std::ops::{Add, Deref, DerefMut};

#[derive(Clone, Debug, PartialEq)]
pub struct BezierCurve<K: Scalar>(pub CurveInternal<K>);
type CurveInternal<K> = SmallVec<[Vector2<K>; 4]>;

impl<K: Scalar> Deref for BezierCurve<K> {
    type Target = CurveInternal<K>;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}
impl<K: Scalar> DerefMut for BezierCurve<K> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl<K: Scalar> BezierCurve<K> {
    pub fn degree(&self) -> usize {
        self.len() - 1
    }
}

impl<K: Num + Scalar> BezierCurve<K> {
    fn usize_to_k(n: usize) -> K {
        let mut k = K::zero();
        for _ in 0..n {
            k = k + K::one();
        }
        k
    }
}

impl<K: ComplexField + Scalar + std::fmt::Display> BezierCurve<K> {
    pub fn elevate(&self) -> BezierCurve<K> {
        let one = K::one();
        let n = BezierCurve::<K>::usize_to_k(self.len());
        let mut points = SmallVec::with_capacity(self.len() + 1);
        let mut j = K::zero(); // Counter converting i from usize to K
        points.push(self[0].clone());
        for i in 1..self.len() {
            j = j.clone() + one.clone();
            let p = &self[i - 1] * (j.clone() / n.clone())
                + &self[i] * ((n.clone() - j.clone()) / n.clone());
            points.push(p);
        }
        points.push(self[self.len() - 1].clone());
        BezierCurve(points)
    }

    pub fn lower(&self) -> BezierCurve<K> {
        let k = self.len();
        let k_as_K = BezierCurve::<K>::usize_to_k(k);
        let n = k - 1;

        assert!(k > 3);

        // build M, which will be (k) rows by (k-1) columns
        let mut M = Matrix::zeros_generic(Dynamic::new(k), Dynamic::new(n));
        M[(0, 0)] = K::one();
        M[(n, n - 1)] = K::one();
        let mut i_as_K = K::zero();
        for i in 1..n {
            i_as_K += K::one();
            M[(i, i - 1)] = i_as_K.clone() / k_as_K.clone();
            M[(i, i)] = K::one() - M[(i, i - 1)].clone();
        }

        // Apply our matrix operations:
        let Mt = M.transpose();
        let mut Mi = M.tr_mul(&M);
        assert!(Mi.try_inverse_mut());
        let V = Mi * Mt;

        // And then we map our k-order list of coordinates
        // to an n-order list of coordinates, instead:
        let old = Matrix2xX::from_columns(&self[..]).transpose();
        let new = &V * &old;
        return BezierCurve(new.row_iter().map(|r| r.transpose()).collect());
    }
}

/* Hitboxes and intersections */
impl<K: RealField + Scalar> BezierCurve<K> {
    pub fn bounding_box(&self) -> BoundingBox<K> {
        BoundingBox::from_slice(&self)
    }

    pub fn minimal_bounding_box(&self) -> BoundingBox<K> {
        assert!(self.degree() < 4);
        let mut points: SmallVec<[Vector2<K>; 6]> = SmallVec::new();
        points.push(self[0].clone());
        points.push(self[self.len() - 1].clone());
        let d = self.derivative();
        let roots_x = Polynomial(d.0.row(0)).roots();
        let roots_y = Polynomial(d.0.row(1)).roots();
        roots_x
            .into_iter()
            .chain(roots_y.into_iter())
            .filter(|t| &K::zero() <= t && t <= &K::one())
            .for_each(|t| {
                points.push(self.castlejau_eval(t));
            });
        BoundingBox::from_iter(points.into_iter())
    }

    pub fn convex_hull(&self) -> Vec<Vector2<K>> {
        convex_hull(self.0.clone().into_vec())
    }

    pub fn locate_point(&self, p: Vector2<K>) -> Option<K> {
        let zero = K::zero();
        let one = K::one();
        let two = one.clone() + one.clone();
        let halve = one.clone() / two.clone();

        #[allow(non_snake_case)]
        let NUM_SUBDIVISIONS: usize = 20;
        #[allow(non_snake_case)]
        let MAX_DEVIATION: K = halve.powi(19);
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
            let f_minus_p: Vector2<K> = function.clone() - p.clone();
            let d_times_d: K = derivative.dot(&derivative);
            let d_times_f_minus_p: K = derivative.dot(&f_minus_p);
            approximation = approximation - d_times_f_minus_p / d_times_d;
            //approximation = approximation - (derivative * (function - p)) / (derivative * derivative);
        }

        Some(approximation)
    }

    pub fn get_intersections(&self, other: &Self) -> Vec<Vector2<K>> {
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
impl<K: Field + Scalar> BezierCurve<K> {
    pub fn split(&self, t: K) -> Option<(BezierCurve<K>, BezierCurve<K>)> {
        /*if t < K::zero() || K::one() < t {
            return None;
        }*/
        if self.len() < 2 {
            return None;
        }
        let inv_t = K::one() - t.clone();
        match &self[..] {
            [a2, b2] => {
                let a1 = a2 * inv_t + b2 * t;
                Some((
                    BezierCurve(smallvec![a2.clone(), a1.clone()]),
                    BezierCurve(smallvec![a1, b2.clone()]),
                ))
            }
            [a3, b3, c3] => {
                let a2 = a3 * inv_t.clone() + b3 * t.clone();
                let b2 = b3 * inv_t.clone() + c3 * t.clone();
                let a1 = &a2 * inv_t + &b2 * t;
                Some((
                    BezierCurve(smallvec![a3.clone(), a2, a1.clone()]),
                    BezierCurve(smallvec![a1, b2, c3.clone()]),
                ))
            }
            [a4, b4, c4, d4] => {
                let a3 = a4 * inv_t.clone() + b4 * t.clone();
                let b3 = b4 * inv_t.clone() + c4 * t.clone();
                let c3 = c4 * inv_t.clone() + d4 * t.clone();
                let a2 = &a3 * inv_t.clone() + &b3 * t.clone();
                let b2 = &b3 * inv_t.clone() + &c3 * t.clone();
                let a1 = &a2 * inv_t + &b2 * t;
                Some((
                    BezierCurve(smallvec![a4.clone(), a3, a2, a1.clone()]),
                    BezierCurve(smallvec![a1, b2, c3, d4.clone()]),
                ))
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
                Some((BezierCurve(lower), BezierCurve(upper)))
            }
        }
    }

    pub fn castlejau_eval(&self, t: K) -> Vector2<K> {
        let inv_t = K::one() - t.clone();
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

    fn castlejau_step(input: &CurveInternal<K>, output: &mut CurveInternal<K>, t: K) {
        output.clear();
        let len = input.len();
        let t_inv = K::one() - t.clone();
        for (p, q) in (&input[0..len - 1]).iter().zip((&input[1..len]).iter()) {
            output.push(p * t_inv.clone() + q * t.clone());
        }
    }
}

/* Stuff using polynomials */
impl<K: Field + Scalar> BezierCurve<K> {
    pub fn polynomial(&self) -> Polynomial2xX<K> {
        let zero = K::zero();
        let one = K::one();
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
                let mut ps = bernstein_polynomials::<K>(self.degree())
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

    pub fn derivative(&self) -> Polynomial2xX<K> {
        let zero = K::zero();
        let one = K::one();
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
                let mut ps = bernstein_polynomials::<K>(self.degree() - 1)
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

    pub fn tangent(&self, t: K) -> Vector2<K> {
        self.derivative().evaluate(t.clone())
    }

    pub fn normal(&self, t: K) -> Vector2<K> {
        let [[x, y]] = self.tangent(t).data.0;
        Vector2::new(K::zero() - y, x)
    }
}
mod convert {
    use nalgebra::{ArrayStorage, Const, Matrix, Matrix2xX, U2};
    pub(crate) fn static_to_dynamic<T: Clone, const C: usize>(
        matrix: &Matrix<T, U2, Const<C>, ArrayStorage<T, 2, C>>,
    ) -> Matrix2xX<T> {
        let mut data = Vec::with_capacity(2 * C);
        data.extend_from_slice(matrix.data.as_slice());
        let data = nalgebra::VecStorage::new(nalgebra::Const::<2>, nalgebra::Dynamic::new(C), data);
        Matrix2xX::from_data(data)
    }
}

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
struct SubCurve<K: RealField> {
    from: K,
    to: K,
    curve: BezierCurve<K>,
}
impl<K: RealField> SubCurve<K> {
    fn split(&self) -> (SubCurve<K>, SubCurve<K>) {
        let two = K::one() + K::one();
        let middle = (self.from.clone() + self.to.clone()) / two.clone();

        let (lower, upper) = self.curve.split(K::one() / two).unwrap();
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

    fn middle_point(&self) -> Vector2<K> {
        let two = K::one() + K::one();
        let middle = (self.from.clone() + self.to.clone()) / two;

        self.curve.castlejau_eval(middle)
    }
}
impl<K: RealField> From<BezierCurve<K>> for SubCurve<K> {
    fn from(curve: BezierCurve<K>) -> Self {
        SubCurve {
            from: K::zero(),
            to: K::one(),
            curve,
        }
    }
}
impl<K: RealField> Deref for SubCurve<K> {
    type Target = BezierCurve<K>;
    fn deref(&self) -> &Self::Target {
        &self.curve
    }
}
impl<K: RealField> DerefMut for SubCurve<K> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.curve
    }
}
