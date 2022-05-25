use std::ops::{Add, AddAssign, Sub, SubAssign, Mul, MulAssign, Div, DivAssign, Deref, DerefMut, Neg};
use num::{Float, Num};

#[derive(Debug, PartialEq, Clone)]
pub struct Polynomial<T: Copy>(pub Vec<T>);

impl <T: Copy> Deref for Polynomial<T> {
    type Target = Vec<T>;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}
impl <T: Copy> DerefMut for Polynomial<T> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl <K: Num + Copy> Polynomial<K> {
    pub fn evaluate(&self, x: K) -> K {
        let mut power = K::one();
        let mut sum = K::zero();
        for &a in self.iter() {
            sum = sum + a * power;
            power = power * x;
        }
        sum
    }

    pub fn derive(&self) -> Polynomial<K> {
        let mut i = K::zero();
        Polynomial(
            self.iter()
                .skip(1)
                .map(|&a| { i = i + K::one(); i * a })
                .collect()
        )
    }
}

impl <K: Float> Polynomial<K> {
    pub fn roots(&self) -> Vec<K> {
        let zero = K::zero();
        let one = K::one();
        let two = one + one;
        let four = two + two;
        match self.as_slice() {
            &[] => Vec::new(),
            &[_] => Vec::new(),
            &[t, m] => vec![(zero - t) / m],
            &[c, b, a] => {
                let d: K = b.powi(2) - four * a * c;
                if d.is_sign_negative() {
                    Vec::new()
                } else if d.is_zero() {
                    vec![(zero - b) / (two * a)]
                } else {
                    let p = zero - b;
                    let q = d.sqrt();
                    let r = two * a;
                    vec![(p + q) / r, (p - q) / r]
                }
            },
            _ => todo!(),
        }
    }
}

macro_rules! impl_componentwise_vector_binary {
    ($Vector:ident, $Trait: ident, $Method:ident, $NegTrait:path, $NegMethod:ident) => {
        impl <'a, K> $Trait<&'a $Vector<K>> for &'a $Vector<K>
            where K: Copy + $Trait<Output=K> + $NegTrait
        {
            type Output = $Vector<K>;
            fn $Method(self, rhs: &'a $Vector<K>) -> Self::Output {
                let mut vec: Vec<K> = self.iter()
                    .zip(rhs.iter())
                    .map(|(&x, &y)| $Trait::$Method(x, y))
                    .collect();

                let len = vec.len();
                for &y in &rhs[len..] {
                    vec.push(y.$NegMethod());
                }
                for &x in &self[len..] {
                    vec.push(x);
                }

                $Vector(vec)
            }
        }
    }
}
impl_componentwise_vector_binary!(Polynomial, Add, add, Clone, clone);
impl_componentwise_vector_binary!(Polynomial, Sub, sub, Neg<Output=K>, neg);

macro_rules! impl_componentwise_vector_assign {
    ($Vector:ident, $Trait: ident, $Method:ident, $NegTrait:path, $NegMethod:ident) => {
        impl <'a, X> $Trait<&'a $Vector<X>> for $Vector<X>
            where X: Copy + $Trait<X> + $NegTrait
        {
            fn $Method(&mut self, rhs: &'a $Vector<X>) {
                self.iter_mut()
                    .zip(rhs.iter())
                    .for_each(|(x, &y)| $Trait::$Method(x, y));
                for y in &rhs[self.len()..] {
                    self.push(y.$NegMethod());
                }
            }
    }
    }
}
impl_componentwise_vector_assign!(Polynomial, AddAssign, add_assign, Clone, clone);
impl_componentwise_vector_assign!(Polynomial, SubAssign, sub_assign, Neg<Output=X>, neg);

impl <'a, X> Mul<&'a Polynomial<X>> for &'a Polynomial<X>
    where X: Copy + Mul<X, Output=X> + Add<X, Output=X>
{
    type Output = Polynomial<X>;
    fn mul(self, rhs: &'a Polynomial<X>) -> Self::Output {
        let mut vec = Vec::with_capacity(self.len() + rhs.len() - 1);

        for (i, &x) in self.iter().enumerate() {
            for (j, &y) in rhs.iter().enumerate() {
                let prod: X = Mul::mul(x, y);
                if let Some(ij) = vec.get_mut(i + j) {
                    *ij = *ij + prod;
                } else {
                    vec.push(prod);
                }
            }
        }

        Polynomial(vec)
    }
}

impl <'a, X> Neg for &'a Polynomial<X>
    where X: Copy + Neg<Output=X>,
{
    type Output = Polynomial<X>;
    fn neg(self) -> Self::Output {
        Polynomial(
            self.iter()
                .map(|&x| Neg::neg(x))
                .collect()
        )
    }
}

macro_rules! impl_componentwise_scalar_binary {
    ($Vector:ident, $Trait: ident, $Method:ident) => {
        impl <'a, X> $Trait<X> for &'a $Vector<X>
            where X: Copy + $Trait<X, Output=X>,
        {
            type Output = $Vector<X>;
            fn $Method(self, rhs: X) -> Self::Output {
                $Vector(
                    self.iter()
                        .map(|&x| $Trait::$Method(x, rhs))
                        .collect()
                )
            }
        }
    }
}
impl_componentwise_scalar_binary!(Polynomial, Mul, mul);
impl_componentwise_scalar_binary!(Polynomial, Div, div);

impl <'a, X> MulAssign<X> for Polynomial<X>
    where X: Copy + Mul<X, Output=X>
{
    fn mul_assign(&mut self, rhs: X) {
        for x in self.iter_mut() {
            *x = *x * rhs;
        }
    }
}

impl <'a, X> DivAssign<X> for Polynomial<X>
    where X: Copy + Div<X, Output=X>
{
    fn div_assign(&mut self, rhs: X) {
        for x in self.iter_mut() {
            *x = *x / rhs;
        }
    }
}
