use std::ops::{Deref, DerefMut, Add};
use nalgebra::{RealField, Vector2, Vector3, Vector4};
use num::{Num, One};
use smallvec::{SmallVec, smallvec};
use crate::bounding_box::BoundingBox;
use crate::polynomial::Polynomial;

#[derive(Clone, Debug, PartialEq)]
pub struct BezierCurve<K: RealField>(pub CurveInternal<K>);
type CurveInternal<K> = SmallVec<[Vector2<K>; 4]>;

impl <K: RealField> Deref for BezierCurve<K> {
    type Target = CurveInternal<K>;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}
impl <K: RealField> DerefMut for BezierCurve<K> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

/* Basic stuff */
impl <K: RealField> BezierCurve<K> {
    pub fn degree(&self) -> usize {
        self.len() - 1
    }

    fn usize_to_k(n: usize) -> K {
        let mut k = K::zero();
        for _ in 0..n {
            k = k + K::one();
        }
        k
    }

    pub fn elevate(&self) -> BezierCurve<K> {
        let one = K::one();
        let n = BezierCurve::<K>::usize_to_k(self.len());
        let mut points = SmallVec::with_capacity(self.len() + 1);
        let mut j = K::zero(); // Counter converting i from usize to K
        points.push(self[0].clone());
        for i in 1..self.len() {
            j = j.clone() + one.clone();
            let p = &self[i - 1] * (j.clone() / n.clone()) + &self[i] * ((n.clone() - j.clone()) / n.clone());
            points.push(p);
        }
        points.push(self[self.len() - 1].clone());
        BezierCurve(points)
    }
}

/* Hitboxes and intersections */
impl <K: RealField> BezierCurve<K> {
    pub fn bounding_box(&self) -> BoundingBox<K> {
        BoundingBox::from_slice(&self)
    }

    pub fn locate_point(&self, p: Vector2<K>) -> Option<K> {
        let zero = K::zero();
        let one = K::one();
        let two = one.clone() + one.clone();
        let halve = one.clone() / two.clone();

        #[allow(non_snake_case)] let NUM_SUBDIVISIONS: usize = 20;
        #[allow(non_snake_case)] let MAX_DEVIATION: K = halve.powi(19);
        #[allow(non_snake_case)] let NUM_NEWTON: usize = 2;

        // Check basic bounding box before doing any allocations
        let bb = self.bounding_box();
        if !bb.contains(p.clone()) {
            return None;
        }

        // Subdivide curve and check bounding boxes
        let mut divisions = (
            &mut vec![SubCurve::from(self.clone())],
            &mut Vec::new(),
        );
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
            deviation = deviation.clone() + (curve.from.clone() - approximation.clone()).powi(2)
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
        #[allow(non_snake_case)] let NUM_SUBDIVISIONS: usize = 20;

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
                    divisions.1.push((s_lower,         o_upper.clone()));
                    divisions.1.push((s_upper.clone(), o_lower));
                    divisions.1.push((s_upper,         o_upper));
                }
            }
            divisions.0.clear();
            divisions = (divisions.1, divisions.0);
        }
        let divisions = divisions.0;
        if divisions.len() == 0 {
            return Vec::new();
        }

        return divisions.iter()
            .map(|(s, _)| s.middle_point())
            .collect()
    }
}

/* Stuff using de castlejau 's algorithm */
impl <K: RealField> BezierCurve<K> {
    pub fn split(&self, t: K) -> Option<(BezierCurve<K>, BezierCurve<K>)> {
        if t < K::zero() || K::one() < t {
            return None;
        }
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
                upper.push(self[len-1].clone());

                let mut points = (
                    &mut self.0.clone(),
                    &mut self.0.clone()
                );
                for i in 1..len {
                    BezierCurve::castlejau_step(points.0, points.1, t.clone());
                    lower.push(points.1[0].clone());
                    upper.push(points.1[len - i - 1].clone());
                    points = (points.1, points.0);
                }
                upper.reverse(); // I find it more intuitive if the t goes through the two parts in the same direction
                Some((
                    BezierCurve(lower),
                    BezierCurve(upper),
                ))
            }
        }
    }

    pub fn castlejau_eval(&self, t: K) -> Vector2<K> {
        let inv_t = K::one() - t.clone();
        match &self[..] {
            [] => panic!(),
            [a1] => a1.clone(),
            [a2, b2] => {
                a2 * inv_t + b2 * t
            },
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

#[derive(Clone)]
struct SubCurve<K: RealField> {
    from: K,
    to: K,
    curve: BezierCurve<K>,
}
impl <K: RealField> SubCurve<K> {
    fn split(&self) -> (SubCurve<K>, SubCurve<K>) {
        let two = K::one() + K::one();
        let middle = (self.from.clone() + self.to.clone()) / two.clone();

        let (lower, upper) = self.curve.split(K::one() / two).unwrap();
        (
            SubCurve {from: self.from.clone(), to: middle.clone(), curve: lower},
            SubCurve {from: middle, to: self.to.clone(), curve: upper}
        )
    }

    fn middle_point(&self) -> Vector2<K> {
        let two = K::one() + K::one();
        let middle = (self.from.clone() + self.to.clone()) / two;

        self.curve.castlejau_eval(middle)
    }
}
impl <K: RealField> From<BezierCurve<K>> for SubCurve<K> {
    fn from(curve: BezierCurve<K>) -> Self {
        SubCurve {from: K::zero(), to: K::one(), curve}
    }
}
impl <K: RealField> Deref for SubCurve<K> {
    type Target = BezierCurve<K>;
    fn deref(&self) -> &Self::Target {
        &self.curve
    }
}
impl <K: RealField> DerefMut for SubCurve<K> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.curve
    }
}

/* Stuff using polynomials */
impl <K: RealField> BezierCurve<K> {
    pub fn x_polynomial(&self) -> Polynomial<K> {
        self.polynomial::<0>()
    }

    pub fn y_polynomial(&self) -> Polynomial<K> {
        self.polynomial::<1>()
    }

    fn polynomial<const I: usize>(&self) -> Polynomial<K> {
        if I > 1 {
            panic!();
        }
        let zero = K::zero();
        let one = K::one();
        match &self[..] {
            [a] => {
                Polynomial(vec![a[I].clone()])
            }
            [a, b] => {
                let p_a = Vector2::new(one.clone(), zero.clone() - one.clone()) * a[I].clone();
                let p_b = Vector2::new(zero, one) * b[I].clone();
                let [p] = (p_a + p_b).data.0;
                Polynomial(p.into())
            }
            [a, b, c] => {
                let two = one.clone() + one.clone();
                let p_a = Vector3::new(one.clone(), zero.clone() - two.clone(), one.clone()) * a[I].clone();
                let p_b = Vector3::new(zero.clone(), two.clone(), zero.clone() - two) * b[I].clone();
                let p_c = Vector3::new(zero.clone(), zero, one) * c[I].clone();
                let [p] = (p_a + p_b + p_c).data.0;
                Polynomial(p.into())
            }
            [a, b, c, d] => {
                let three = one.clone() + one.clone() + one.clone();
                let six = three.clone() + three.clone();
                let p_a = Vector4::new(one.clone(), zero.clone() - three.clone(), three.clone(), zero.clone() - one.clone()) * a[I].clone();
                let p_b = Vector4::new(zero.clone(), three.clone(), zero.clone() - six, three.clone()) * b[I].clone();
                let p_c = Vector4::new(zero.clone(), zero.clone(), three.clone(), zero.clone() - three) * c[I].clone();
                let p_d = Vector4::new(zero.clone(), zero.clone(), zero, one) * d[I].clone();
                let [p] = (p_a + p_b + p_c + p_d).data.0;
                Polynomial(p.into())
            }
            _ => {
                let mut ps = bernstein_polynomials::<K>(self.degree())
                    .into_iter()
                    .zip(self.iter())
                    .map(|(mut p, vec)| {
                        let x = &vec[I];
                        for a in p.iter_mut() {
                            *a = a.clone() * x.clone();
                        }
                        p
                    });
                if let Some(mut p) = ps.next() {
                    for q in ps {
                        // Manually pasted and adjusted AddAssign
                        p.iter_mut()
                            .zip(q.iter())
                            .for_each(|(x, y)| *x = x.clone() + y.clone());
                        for y in &q[p.len()..] {
                            p.push(y.clone());
                        }
                    }
                    p
                } else {
                    unreachable!();
                }
            }
        }
    }

    pub fn x_derivative(&self) -> Polynomial<K> {
        self.derivative::<0>()
    }

    pub fn y_derivative(&self) -> Polynomial<K> {
        self.derivative::<1>()
    }

    fn derivative<const I: usize>(&self) -> Polynomial<K> {
        if I > 1 {
            panic!();
        }
        let zero = K::zero();
        let one = K::one();
        match &self[..] {
            [_] => {
                Polynomial(vec![])
            }
            [a, b] => {
                Polynomial(vec![(b - a)[I].clone()])
            }
            [a, b, c] => {
                let two = one.clone() + one.clone();
                let p_a = Vector2::new(one.clone(), zero.clone() - one.clone()) * (b - a)[I].clone();
                let p_b = Vector2::new(zero, one) * (c - b)[I].clone();
                let [p] = ((p_a + p_b) * two).data.0;
                Polynomial(p.into())
            },
            [a, b, c, d] => {
                let two = one.clone() + one.clone();
                let three = two.clone() + one.clone();
                let p_a = Vector3::new(one.clone(), zero.clone() - two.clone(), one.clone()) * (b - a)[I].clone();
                let p_b = Vector3::new(zero.clone(), two.clone(), zero.clone() - two) * (c - b)[I].clone();
                let p_c = Vector3::new(zero.clone(), zero, one) * (d - c)[I].clone();
                let [p] = ((p_a + p_b + p_c) * three).data.0;
                Polynomial(p.into())
            }
            _ => {
                let mut degree = zero;
                let mut ps = bernstein_polynomials::<K>(self.degree()-1)
                    .into_iter()
                    .enumerate()
                    .map(|(i, mut p)| {
                        degree = degree.clone() + one.clone();
                        let point = &self[i+1] - &self[i];
                        for a in p.iter_mut() {
                            *a = a.clone() * point[I].clone();
                        }
                        p
                    });
                if let Some(mut p) = ps.next() {
                    for q in ps {
                        // Manually pasted and adjusted AddAssign
                        p.iter_mut()
                            .zip(q.iter())
                            .for_each(|(x, y)| *x = x.clone() + y.clone());
                        for y in &q[p.len()..] {
                            p.push(y.clone());
                        }
                    }
                    &p * degree
                } else {
                    unreachable!();
                }
            }
        }
    }

    pub fn tangent(&self, t: K) -> Vector2<K> {
        Vector2::new(
            self.x_derivative().evaluate(t.clone()),
            self.y_derivative().evaluate(t),
        )
    }

    pub fn normal(&self, t: K) -> Vector2<K> {
        let [[x, y]] = self.tangent(t).data.0;
        Vector2::new(
            K::zero() - y,
            x,
        )
    }

    pub fn minimal_bounding_box(&self) -> BoundingBox<K> {
        assert!(self.degree() < 4);
        let mut points: SmallVec<[Vector2<K>; 6]> = SmallVec::new();
        points.push(self[0].clone());
        points.push(self[self.len()-1].clone());
        self.x_derivative().roots().into_iter()
            .chain(self.y_derivative().roots().into_iter())
            .filter(|t| &K::zero() <= t && t <= &K::one())
            .for_each(|t| {
                points.push(self.castlejau_eval(t));
            });
        BoundingBox::from_iter(points.into_iter())
    }
}

pub fn bernstein_polynomials<N>(degree: usize) -> Vec<Polynomial<N>>
    where N: Num + Clone
{
    let one = N::one();
    let zero = N::zero();

    // Prepare the powers of x and (1-x)
    let mut powers = (
        Vec::with_capacity(degree+1),
        Vec::with_capacity(degree+1),
    );

    powers.0.push(Polynomial(vec![one.clone()]));
    powers.1.push(Polynomial(vec![one.clone()]));
    powers.0.push(Polynomial(vec![zero.clone(), one.clone()]));
    powers.1.push(Polynomial(vec![one.clone(), zero - one]));

    for i in 1..degree {
        powers.0.push(&powers.0[i] * &powers.0[1]);
        powers.1.push(&powers.1[i] * &powers.1[1]);
    }

    // Combine powers into Bernstein polynomials
    let mut base = Vec::with_capacity(degree+1);
    let pascal = pascal_triangle::<N>(degree);
    for i in 0..degree+1 {
        let mut b = &powers.0[i] * &powers.1[degree-i];
        b *= pascal[i].clone(); // TODO use into_iter
        base.push(b);
    }

    return base;
}

pub fn pascal_triangle<N>(layer: usize) -> Vec<N>
    where N: Add<Output=N> + One + Clone
{
    let one = N::one();
    let mut old_layer = Vec::with_capacity(layer+1);
    let mut new_layer = Vec::with_capacity(layer+1);
    new_layer.push(N::one());

    while new_layer.len() < new_layer.capacity() {
        old_layer = new_layer.clone(); // TODO fix performance

        new_layer.push(one.clone());
        let get = |i| old_layer.get(i as usize).map(|n| n).unwrap_or(&one).clone();
        for i in 1..new_layer.len()-1 {
            new_layer[i] = get(i-1) + get(i);
        }
    }

    new_layer
}