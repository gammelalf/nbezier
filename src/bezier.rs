use std::ops::{Deref, DerefMut, AddAssign};
use num::Float;
use crate::polynomial::Polynomial;
use crate::vector::Vector;

impl <'a, K: Float, const N: usize> AddAssign<&'a Vector<K, N>> for Vector<K, N> {
    fn add_assign(&mut self, rhs: &Vector<K, N>) {
        for (x, &y) in self.iter_mut().zip(rhs.iter()) {
            *x = *x + y;
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct BezierCurve<K: Float>(pub Vec<Vector<K, 2>>);
impl <K: Float> Deref for BezierCurve<K> {
    type Target = Vec<Vector<K, 2>>;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}
impl <K: Float> DerefMut for BezierCurve<K> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

/* Stuff stuff */
impl <K: Float> BezierCurve<K> {
    pub fn bounding_box(&self) -> [Vector<K, 2>; 2] {
        let mut min = self[0];
        let mut max = self[0];
        for p in self.iter().skip(1) {
            if min[0] > p[0] {
                min[0] = p[0];
            }
            if min[1] > p[1] {
                min[1] = p[1];
            }
            if max[0] < p[0] {
                max[0] = p [0];
            }
            if max[1] < p[1] {
                max[1] = p[1];
            }
        }
        return [min, max];
    }
}

/* Stuff using de castlejau 's algorithm */
impl <K: Float> BezierCurve<K> {
    pub fn degree(&self) -> usize {
        self.len() - 1
    }

    pub fn split(&self, t: K) -> Option<(BezierCurve<K>, BezierCurve<K>)> {
        if t < K::zero() || K::one() < t {
            None
        } else if self.len() < 2 {
            None
        } else {
            let len = self.len();

            let mut lower = Vec::with_capacity(len);
            let mut upper = Vec::with_capacity(len);
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

    pub fn castlejau_eval(&self, t: K) -> Vector<K, 2> {
        let mut old_points = self.0.clone();
        let mut new_points = self.0.clone();
        let mut points = (&mut old_points, &mut new_points);
        while points.1.len() > 1 {
            BezierCurve::castlejau_step(points.0, points.1, t);
            points = (points.1, points.0);
        }
        return points.1[0];
    }

    fn castlejau_step(input: &Vec<Vector<K, 2>>, output: &mut Vec<Vector<K, 2>>, t: K) {
        output.clear();
        let len = input.len();
        let t_inv = K::one() - t;
        for (&p, &q) in (&input[0..len - 1]).iter().zip((&input[1..len]).iter()) {
            output.push(p * t_inv + q * t);
        }
    }
}

/* Stuff using polynomials */
impl <K: Float + From<i32>> BezierCurve<K> {
    pub fn x_polynomial(&self) -> Polynomial<K> {
        let ps = get_bernstein_polynomials(self.degree())
            .into_iter()
            .zip(self.iter())
            .map(|(p, Vector([x, _]))| {
                let mut p: Polynomial<K> = p.convert();
                for a in p.iter_mut() {
                    *a = *a * *x;
                }
                p
            });
        BezierCurve::sum_polynomials(ps)
    }

    pub fn y_polynomial(&self) -> Polynomial<K> {
        let ps = get_bernstein_polynomials(self.degree())
            .into_iter()
            .zip(self.iter())
            .map(|(p, Vector([_, y]))| {
                let mut p: Polynomial<K> = p.convert();
                for a in p.iter_mut() {
                    *a = *a * *y;
                }
                p
            });
        BezierCurve::sum_polynomials(ps)
    }

    fn sum_polynomials<I: Iterator<Item=Polynomial<K>>>(ps: I) -> Polynomial<K> {
        let mut ps = ps;
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
            Polynomial(vec![])
        }
    }

    pub fn tangent(&self, t: K) -> Vector<K, 2> {
        Vector([
            self.x_polynomial().derive().evaluate(t),
            self.y_polynomial().derive().evaluate(t),
        ])
    }

    pub fn normal(&self, t: K) -> Vector<K, 2> {
        let tangent = self.tangent(t);
        Vector([
            K::zero() - tangent[1],
            tangent[0],
        ])
    }
}

pub fn get_bernstein_polynomials(degree: usize) -> Vec<Polynomial<i32>> {
    #[cfg(feature = "cache")]
    if degree < 4 {
        return BERNSTEIN_POLYNOMIALS[degree].clone();
    }

    return calc_bernstein_polynomials(degree);
}

pub fn calc_bernstein_polynomials(degree: usize) -> Vec<Polynomial<i32>> {
    // Prepare the powers of x and (1-x)
    let mut powers = (
        Vec::with_capacity(degree+1),
        Vec::with_capacity(degree+1),
    );

    powers.0.push(Polynomial(vec![1]));
    powers.1.push(Polynomial(vec![1]));
    powers.0.push(Polynomial(vec![0,   1]));
    powers.1.push(Polynomial(vec![1,  -1]));

    for i in 1..degree {
        powers.0.push(&powers.0[i] * &powers.0[1]);
        powers.1.push(&powers.1[i] * &powers.1[1]);
    }

    // Combine powers into Bernstein polynomials
    let mut base = Vec::with_capacity(degree+1);
    let pascal = pascal_triangle(degree);
    for i in 0..degree+1 {
        let mut b = &powers.0[i] * &powers.1[degree-i];
        b *= pascal[i];
        base.push(b);
    }

    return base;
}

pub fn pascal_triangle(layer: usize) -> Vec<i32> {
    let mut old_layer = Vec::with_capacity(layer+1);
    let mut new_layer = Vec::with_capacity(layer+1);
    new_layer.push(1);

    while new_layer.len() < new_layer.capacity() {
        old_layer.push(0);
        old_layer.copy_from_slice(new_layer.as_slice());

        new_layer.push(1);
        let get = |i| old_layer.get(i as usize).map(|&n| n).unwrap_or(1);
        for i in 1..new_layer.len()-1 {
            new_layer[i] = get(i-1) + get(i);
        }
    }

    new_layer
}

#[cfg(feature = "cache")]
use once_cell::sync::Lazy;

#[cfg(feature = "cache")]
static BERNSTEIN_POLYNOMIALS: [Lazy<Vec<Polynomial<i32>>>; 4] = [
    Lazy::new(|| calc_bernstein_polynomials(0)),
    Lazy::new(|| calc_bernstein_polynomials(1)),
    Lazy::new(|| calc_bernstein_polynomials(2)),
    Lazy::new(|| calc_bernstein_polynomials(3)),
];