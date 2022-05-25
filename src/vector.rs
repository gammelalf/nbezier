use std::cmp::Ordering;
use std::ops::{Add, Sub, Mul, Div, Deref, DerefMut};
use num::{Float, Zero};

#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Vector<K: Copy, const N: usize>(pub [K; N]);

/* Deref */
impl <K: Copy, const N: usize> Deref for Vector<K, N> {
    type Target = [K; N];
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}
impl <K: Copy, const N: usize> DerefMut for Vector<K, N> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

/* Compare points suitable to for aabbs */
impl <K: Copy, const N: usize> PartialOrd for Vector<K, N> where K: PartialOrd {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        let mut ordering = None;
        for (x, y) in self.iter().zip(other.iter()) {
            match x.partial_cmp(y) {
                None => return None, // If a coordinate can't be compared, the whole vector can't be
                Some(Ordering::Equal) => (), // Equal coordinates are just skipped
                Some(next_ordering) => {
                    if let Some(old_ordering) = ordering {
                        if next_ordering != old_ordering {
                            return None; // A coordinate is less while another is greater => can't be compared
                        }
                    } else {
                        ordering = Some(next_ordering);
                    }
                }
            }
        }
        if ordering.is_none() {
            Some(Ordering::Equal)
        } else {
            ordering
        }
    }
}

/* Basic arithmetic */
impl <K: Copy, const N: usize> Add for Vector<K, N> where K: Add<K, Output=K> {
    type Output = Vector<K, N>;
    fn add(mut self, rhs: Self) -> Self::Output {
        self.iter_mut()
            .zip(rhs.into_iter())
            .for_each(|(x, y)| *x = Add::add(*x, y));
        self
    }
}
impl<K: Copy, const N: usize> Sub for Vector<K, N> where K: Sub<K, Output=K> {
    type Output = Vector<K, N>;
    fn sub(mut self, rhs: Self) -> Self::Output {
        self.iter_mut()
            .zip(rhs.into_iter())
            .for_each(|(x, y)| *x = Sub::sub(*x, y));
        self
    }
}

impl <K: Copy, const N: usize> Mul<K> for Vector<K, N> where K: Mul<K, Output=K> {
    type Output = Vector<K, N>;
    fn mul(mut self, y: K) -> Self::Output {
        self.iter_mut()
            .for_each(|x| *x = Mul::mul(*x, y));
        self
    }
}
impl <K: Copy, const N: usize> Div<K> for Vector<K, N> where K: Div<K, Output=K> {
    type Output = Vector<K, N>;
    fn div(mut self, y: K) -> Self::Output {
        self.iter_mut()
            .for_each(|x| *x = Div::div(*x, y));
        self
    }
}

/* Scalar product and euclidean norm */
impl <K: Copy, const N: usize> Mul for Vector<K, N> where K: Mul<K, Output=K> + Add<K, Output=K> + Zero {
    type Output = K;
    fn mul(self, rhs: Self) -> Self::Output {
        self.into_iter()
            .zip(rhs.into_iter())
            .map(|(x, y)| x * y)
            .fold(K::zero(), |x, y| x + y)
    }
}
impl <K: Copy, const N: usize> Vector<K, N> where K: Float {
    pub fn norm(&self) -> K {
        (*self * *self).sqrt()
    }

    pub fn normalize(&self) -> Self {
        *self / self.norm()
    }

    pub fn angle(&self, &other: Self) -> K {
        ((*self * *other) / (self.norm() * other.norm())).acos()
    }
}

/* Cross product */
impl <K: Copy> Vector<K, 2> where K: Mul<K, Output=K> + Sub<K, Output=K> {
    /* Compute the would be third component */
    pub fn cross(&self, &other: Self) -> K {
        self[0] * other[1] - self[1] * other[0]
    }
}
impl <K: Copy> Vector<K, 3> where K: Mul<K, Output=K> + Sub<K, Output=K> {
    pub fn cross(&self, &other: Self) -> Vector<K, 3> {
        Vector([
            self[1] * other[2] - self[2] * other[1],
            self[2] * other[0] - self[0] * other[2],
            self[0] * other[1] - self[1] * other[0],
        ])
    }
}

/* Conversions from and to tuples */
impl <K: Copy> From<(K, K)> for Vector<K, 2> {
    fn from(p: (K, K)) -> Self {
        Vector([p.0, p.1])
    }
}
impl <K: Copy> From<Vector<K, 2>> for (K, K) {
    fn from(p: Vector<K, 2>) -> Self {
        (p[0], p[1])
    }
}
impl <K: Copy> From<(K, K, K)> for Vector<K, 3> {
    fn from(p: (K, K, K)) -> Self {
        Vector([p.0, p.1, p.2])
    }
}
impl <K: Copy> From<Vector<K, 3>> for (K, K, K) {
    fn from(p: Vector<K, 3>) -> Self {
        (p[0], p[1], p[2])
    }
}
