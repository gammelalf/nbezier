use std::ops::{Deref, DerefMut, Add};
use num::{Float, Num, One};
use smallvec::{SmallVec, smallvec};
use crate::vector::Vector;
use crate::bounding_box::BoundingBox;
use crate::polynomial::Polynomial;

#[derive(Clone, Debug, PartialEq)]
pub struct BezierCurve<K: Float>(pub CurveInternal<K>);
type CurveInternal<K> = SmallVec<[Vector<K, 2>; 4]>;

impl <K: Float> Deref for BezierCurve<K> {
    type Target = CurveInternal<K>;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}
impl <K: Float> DerefMut for BezierCurve<K> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

/* Basic stuff */
impl <K: Float> BezierCurve<K> {
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
        points.push(self[0]);
        for i in 1..self.len() {
            j = j + one;
            let p = self[i - 1] * j / n + self[i] * (n - j) / n;
            points.push(p);
        }
        points.push(self[self.len() - 1]);
        BezierCurve(points)
    }
}

/* Hitboxes and intersections */
impl <K: Float> BezierCurve<K> {
    pub fn bounding_box(&self) -> BoundingBox<K> {
        BoundingBox::from_slice(&self)
    }

    pub fn locate_point(&self, p: Vector<K, 2>) -> Option<K> {
        let zero = K::zero();
        let one = K::one();
        let two = one + one;
        let halve = one / two;

        #[allow(non_snake_case)] let NUM_SUBDIVISIONS: usize = 20;
        #[allow(non_snake_case)] let MAX_DEVIATION: K = halve.powi(19);
        #[allow(non_snake_case)] let NUM_NEWTON: usize = 2;

        // Check basic bounding box before doing any allocations
        let bb = self.bounding_box();
        if !bb.contains(p) {
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
                if bb.contains(p) {
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
        let mut divisions_len = zero; // Replacement for divisions.len() which would be usize
        let mut approximation = zero;
        for curve in divisions.iter() {
            approximation = approximation + curve.from + curve.to;
            divisions_len = divisions_len + one;
        }
        approximation = approximation / (two * divisions_len);

        // Check deviation to see if subdivisions are too far apart
        let mut deviation = zero;
        for curve in divisions.iter() {
            deviation = deviation + (curve.from - approximation).powi(2)
                                  + (curve.to - approximation).powi(2);
        }
        deviation = (deviation / (two * divisions_len)).sqrt();
        if deviation > MAX_DEVIATION {
            return None;
        }

        // Refine solution using newton's method
        for _ in 0..NUM_NEWTON {
            let function = self.castlejau_eval(approximation);
            let derivative = self.tangent(approximation);
            approximation = approximation - (derivative * (function - p)) / (derivative * derivative);
        }

        Some(approximation)
    }

    pub fn get_intersections(&self, other: &Self) -> Vec<Vector<K, 2>> {
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
impl <K: Float> BezierCurve<K> {
    pub fn split(&self, t: K) -> Option<(BezierCurve<K>, BezierCurve<K>)> {
        if t < K::zero() || K::one() < t {
            return None;
        }
        if self.len() < 2 {
            return None;
        }
        let inv_t = K::one() - t;
        match &self[..] {
            &[a2, b2] => {
                let a1 = a2 * inv_t + b2 * t;
                Some((
                    BezierCurve(smallvec![a2, a1]),
                    BezierCurve(smallvec![a1, b2]),
                ))
            }
            &[a3, b3, c3] => {
                let a2 = a3 * inv_t + b3 * t;
                let b2 = b3 * inv_t + c3 * t;
                let a1 = a2 * inv_t + b2 * t;
                Some((
                    BezierCurve(smallvec![a3, a2, a1]),
                    BezierCurve(smallvec![a1, b2, c3]),
                ))
            }
            &[a4, b4, c4, d4] => {
                let a3 = a4 * inv_t + b4 * t;
                let b3 = b4 * inv_t + c4 * t;
                let c3 = c4 * inv_t + d4 * t;
                let a2 = a3 * inv_t + b3 * t;
                let b2 = b3 * inv_t + c3 * t;
                let a1 = a2 * inv_t + b2 * t;
                Some((
                    BezierCurve(smallvec![a4, a3, a2, a1]),
                    BezierCurve(smallvec![a1, b2, c3, d4]),
                ))
            }
            _ => {
                let len = self.len();

                let mut lower = SmallVec::with_capacity(len);
                let mut upper = SmallVec::with_capacity(len);
                lower.push(self[0]);
                upper.push(self[len-1]);

                let mut old_points = self.0.clone();
                let mut new_points = self.0.clone();
                let mut points = (&mut old_points, &mut new_points);
                for i in 1..len {
                    BezierCurve::castlejau_step(points.0, points.1, t);
                    lower.push(points.1[0]);
                    upper.push(points.1[len - i - 1]);
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

    pub fn castlejau_eval(&self, t: K) -> Vector<K, 2> {
        let inv_t = K::one() - t;
        match &self[..] {
            &[] => panic!(),
            &[a1] => a1,
            &[a2, b2] => {
                a2 * inv_t + b2 * t
            },
            &[a3, b3, c3] => {
                let a2 = a3 * inv_t + b3 * t;
                let b2 = b3 * inv_t + c3 * t;
                a2 * inv_t + b2 * t
            }
            &[a4, b4, c4, d4] => {
                let a3 = a4 * inv_t + b4 * t;
                let b3 = b4 * inv_t + c4 * t;
                let c3 = c4 * inv_t + d4 * t;
                let a2 = a3 * inv_t + b3 * t;
                let b2 = b3 * inv_t + c3 * t;
                a2 * inv_t + b2 * t
            }
            _ => {
                let mut old_points = self.0.clone();
                let mut new_points = self.0.clone();
                let mut points = (&mut old_points, &mut new_points);
                while points.1.len() > 1 {
                    BezierCurve::castlejau_step(points.0, points.1, t);
                    points = (points.1, points.0);
                }
                return points.1[0];
            }
        }
    }

    fn castlejau_step(input: &CurveInternal<K>, output: &mut CurveInternal<K>, t: K) {
        output.clear();
        let len = input.len();
        let t_inv = K::one() - t;
        for (&p, &q) in (&input[0..len - 1]).iter().zip((&input[1..len]).iter()) {
            output.push(p * t_inv + q * t);
        }
    }
}

#[derive(Clone)]
struct SubCurve<K: Float> {
    from: K,
    to: K,
    curve: BezierCurve<K>,
}
impl <K: Float> SubCurve<K> {
    fn split(&self) -> (SubCurve<K>, SubCurve<K>) {
        let two = K::one() + K::one();
        let middle = (self.from + self.to) / two;

        let (lower, upper) = self.curve.split(K::one() / two).unwrap();
        (
            SubCurve {from: self.from, to: middle, curve: lower},
            SubCurve {from: middle, to: self.to, curve: upper}
        )
    }

    fn middle_point(&self) -> Vector<K, 2> {
        let two = K::one() + K::one();
        let middle = (self.from + self.to) / two;

        self.curve.castlejau_eval(middle)
    }
}
impl <K: Float> From<BezierCurve<K>> for SubCurve<K> {
    fn from(curve: BezierCurve<K>) -> Self {
        SubCurve {from: K::zero(), to: K::one(), curve}
    }
}
impl <K: Float> Deref for SubCurve<K> {
    type Target = BezierCurve<K>;
    fn deref(&self) -> &Self::Target {
        &self.curve
    }
}
impl <K: Float> DerefMut for SubCurve<K> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.curve
    }
}

/* Stuff using polynomials */
impl <K: Float> BezierCurve<K> {
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
            &[a] => {
                Polynomial(vec![a[I]])
            }
            &[a, b] => {
                let p_a = Vector([one, zero - one]) * a[I];
                let p_b = Vector([zero, one]) * b[I];
                Polynomial(vec![
                    p_a[0] + p_b[0],
                    p_a[1] + p_b[1],
                ])
            }
            &[a, b, c] => {
                let two = one + one;
                let p_a = Vector([one, zero - two, one]) * a[I];
                let p_b = Vector([zero, two, zero - two]) * b[I];
                let p_c = Vector([zero, zero, one]) * c[I];
                Polynomial(vec![
                    p_a[0] + p_b[0] + p_c[0],
                    p_a[1] + p_b[1] + p_c[1],
                    p_a[2] + p_b[2] + p_c[2],
                ])
            }
            &[a, b, c, d] => {
                let three = one + one + one;
                let six = three + three;
                let p_a = Vector([one, zero - three, three, zero - one]) * a[I];
                let p_b = Vector([zero, three, zero - six, three]) * b[I];
                let p_c = Vector([zero, zero, three, zero - three]) * c[I];
                let p_d = Vector([zero, zero, zero, one]) * d[I];
                Polynomial(vec![
                    p_a[0] + p_b[0] + p_c[0] + p_d[0],
                    p_a[1] + p_b[1] + p_c[1] + p_d[1],
                    p_a[2] + p_b[2] + p_c[2] + p_d[2],
                    p_a[3] + p_b[3] + p_c[3] + p_d[3],
                ])
            }
            _ => {
                let mut ps = bernstein_polynomials::<K>(self.degree())
                    .into_iter()
                    .zip(self.iter())
                    .map(|(mut p, Vector([x, _]))| {
                        for a in p.iter_mut() {
                            *a = *a * *x;
                        }
                        p
                    });
                if let Some(mut p) = ps.next() {
                    for q in ps {
                        // Manually pasted and adjusted AddAssign
                        p.iter_mut()
                            .zip(q.iter())
                            .for_each(|(x, y)| *x = *x + *y);
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
            &[_] => {
                Polynomial(vec![])
            }
            &[a, b] => {
                Polynomial(vec![(b - a)[I]])
            }
            &[a, b, c] => {
                let two = one + one;
                let p_a = Vector([one, zero - one]) * (b - a)[I];
                let p_b = Vector([zero, one]) * (c - b)[I];
                Polynomial(vec![
                    two * (p_a[0] + p_b[0]),
                    two * (p_a[1] + p_b[1]),
                ])
            },
            &[a, b, c, d] => {
                let two = one + one;
                let three = two + one;
                let p_a = Vector([one, zero - two, one]) * (b - a)[I];
                let p_b = Vector([zero, two, zero - two]) * (c - b)[I];
                let p_c = Vector([zero, zero, one]) * (d - c)[I];
                Polynomial(vec![
                    three * (p_a[0] + p_b[0] + p_c[0]),
                    three * (p_a[1] + p_b[1] + p_c[1]),
                    three * (p_a[2] + p_b[2] + p_c[2]),
                ])
            }
            _ => {
                let mut degree = zero;
                let mut ps = bernstein_polynomials::<K>(self.degree()-1)
                    .into_iter()
                    .enumerate()
                    .map(|(i, mut p)| {
                        degree = degree + one;
                        let point = self[i+1] - self[i];
                        for a in p.iter_mut() {
                            *a = *a * point[I];
                        }
                        p
                    });
                if let Some(mut p) = ps.next() {
                    for q in ps {
                        // Manually pasted and adjusted AddAssign
                        p.iter_mut()
                            .zip(q.iter())
                            .for_each(|(x, y)| *x = *x + *y);
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

    pub fn tangent(&self, t: K) -> Vector<K, 2> {
        Vector([
            self.x_derivative().evaluate(t),
            self.y_derivative().evaluate(t),
        ])
    }

    pub fn normal(&self, t: K) -> Vector<K, 2> {
        let tangent = self.tangent(t);
        Vector([
            K::zero() - tangent[1],
            tangent[0],
        ])
    }

    pub fn minimal_bounding_box(&self) -> BoundingBox<K> {
        assert!(self.degree() < 4);
        let mut points: SmallVec<[Vector<K, 2>; 6]> = SmallVec::new();
        points.push(self[0]);
        points.push(self[self.len()-1]);
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
    where N: Num + Copy
{
    let one = N::one();
    let zero = N::zero();

    // Prepare the powers of x and (1-x)
    let mut powers = (
        Vec::with_capacity(degree+1),
        Vec::with_capacity(degree+1),
    );

    powers.0.push(Polynomial(vec![one]));
    powers.1.push(Polynomial(vec![one]));
    powers.0.push(Polynomial(vec![zero, one]));
    powers.1.push(Polynomial(vec![one, zero - one]));

    for i in 1..degree {
        powers.0.push(&powers.0[i] * &powers.0[1]);
        powers.1.push(&powers.1[i] * &powers.1[1]);
    }

    // Combine powers into Bernstein polynomials
    let mut base = Vec::with_capacity(degree+1);
    let pascal = pascal_triangle::<N>(degree);
    for i in 0..degree+1 {
        let mut b = &powers.0[i] * &powers.1[degree-i];
        b *= pascal[i];
        base.push(b);
    }

    return base;
}

pub fn pascal_triangle<N>(layer: usize) -> Vec<N>
    where N: Add<Output=N> + One + Copy
{
    let one = N::one();
    let mut old_layer = Vec::with_capacity(layer+1);
    let mut new_layer = Vec::with_capacity(layer+1);
    new_layer.push(N::one());

    while new_layer.len() < new_layer.capacity() {
        old_layer.push(one);
        old_layer.copy_from_slice(new_layer.as_slice());

        new_layer.push(one);
        let get = |i| old_layer.get(i as usize).map(|&n| n).unwrap_or(one);
        for i in 1..new_layer.len()-1 {
            new_layer[i] = get(i-1) + get(i);
        }
    }

    new_layer
}