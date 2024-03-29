//! A wrapper around [`nalgebra::Matrix`] interpreting it as a polynomial.

use nalgebra::allocator::Allocator;
use nalgebra::dimension::{Const, Dim, DimAdd, DimDiff, DimSub, DimSum, U1, U2, U3, U4};
use nalgebra::storage::{RawStorage, Storage, StorageMut};
use nalgebra::{
    DefaultAllocator, Dynamic, Field, Matrix, OMatrix, OVector, Owned, RealField, Scalar,
};
use std::fmt::{self, Write};

/// Wrapper around [`nalgebra::OMatrix`] interpreting it as a polynomial.
pub type OPolynomial<T, R, C> = Polynomial<T, R, C, Owned<T, R, C>>;

/// Wrapper around [`nalgebra::Matrix`] interpreting it as a polynomial:
/// $p: \R \to \R^r $ where $r$ is the number of rows i.e. the generic `R` parameter
///
/// This means rows are the polynomials for each coordinate
/// and columns are the different powers' coefficents.
pub struct Polynomial<T, R, C, S>(pub Matrix<T, R, C, S>);

impl<T: RealField, R: Dim, C: Dim, S: Storage<T, R, C>> Polynomial<T, R, C, S> {
    /// Get the polynomial degree i.e its highest power
    ///
    /// For example a quadratic polynomial has degree 2 and 3 coefficients
    pub fn degree(&self) -> usize {
        self.0.ncols() - 1
    }
}

/* Eval, Derive, Integrate */
impl<T: Scalar, R: Dim, C: Dim, S: Storage<T, R, C>> Polynomial<T, R, C, S> {
    /// Evaluate `self` at position `x` and store the result into `out`.
    pub fn evaluate_to<S2>(&self, x: T, out: &mut Matrix<T, R, U1, S2>)
    where
        T: Field,
        S2: StorageMut<T, R, U1>,
    {
        out.fill(T::zero());
        for i in 0..self.0.ncols() - 1 {
            let i = self.0.ncols() - 1 - i;
            *out += self.0.column(i);
            *out *= x.clone();
        }
        *out += self.0.column(0);
    }

    /// Evaluate `self` at position `x`.
    pub fn evaluate(&self, x: T) -> OVector<T, R>
    where
        T: Field,
        DefaultAllocator: Allocator<T, R>,
    {
        let (rows, _) = self.0.shape_generic();
        let mut out = OVector::from_element_generic(rows, Const::<1>, T::zero());
        self.evaluate_to(x, &mut out);
        out
    }

    /// Calculate `self`'s derivative and store the result into `out`.
    pub fn derive_to<S2>(&self, out: &mut Matrix<T, R, DimDiff<C, U1>, S2>)
    where
        T: Field,
        C: DimSub<U1>,
        S2: StorageMut<T, R, DimDiff<C, U1>>,
    {
        let mut exponent = T::one();
        for i in 1..self.0.ncols() {
            out.set_column(i - 1, &self.0.column(i));
            *&mut out.column_mut(i - 1) *= exponent.clone();
            exponent += T::one();
        }
    }

    /// Calculate `self`'s derivative.
    pub fn derive(&self) -> Polynomial<T, R, DimDiff<C, U1>, Owned<T, R, DimDiff<C, U1>>>
    where
        T: Field,
        C: DimSub<U1>,
        DefaultAllocator: Allocator<T, R, DimDiff<C, U1>>,
    {
        let (c, r) = self.0.shape_generic();
        let mut out = OMatrix::zeros_generic(c, r.sub(Const::<1>));
        self.derive_to(&mut out);
        Polynomial(out)
    }

    /// Calculate `self`'s integral and store the result into `out`.
    ///
    /// `out`'s first column, i.e. what would be the integration constant, is left untouched.
    pub fn integrate_to<S2>(&self, out: &mut Matrix<T, R, DimSum<C, U1>, S2>)
    where
        T: Field,
        C: DimAdd<U1>,
        S2: StorageMut<T, R, DimSum<C, U1>>,
    {
        let mut exponent = T::one();
        for i in 0..self.0.ncols() {
            out.set_column(i + 1, &self.0.column(i));
            *&mut out.column_mut(i + 1) /= exponent.clone();
            exponent += T::one();
        }
    }

    /// Calculate `self`'s integral.
    ///
    /// The zero vector is used as integration constant.
    pub fn integrate(&self) -> Polynomial<T, R, DimSum<C, U1>, Owned<T, R, DimSum<C, U1>>>
    where
        T: Field,
        C: DimAdd<U1>,
        DefaultAllocator: Allocator<T, R, DimSum<C, U1>>,
    {
        let (r, c) = self.0.shape_generic();
        let mut out = OMatrix::zeros_generic(r, c.add(Const::<1>));
        self.integrate_to(&mut out);
        Polynomial(out)
    }
}

/* Product */
impl<T: Scalar, R: Dim, C: Dim, S: Storage<T, R, C>> Polynomial<T, R, C, S> {
    /// Mutliply `self` with `rhs`.
    pub fn mul<CR, SR>(
        &self,
        rhs: &Matrix<T, R, CR, SR>,
    ) -> Polynomial<T, R, DimPolyProd<C, CR>, Owned<T, R, DimPolyProd<C, CR>>>
    where
        T: Field,
        CR: Dim,
        C: DimAdd<CR>,
        DimSum<C, CR>: DimSub<U1>,
        SR: Storage<T, R, CR>,
        DefaultAllocator: Allocator<T, R, DimPolyProd<C, CR>>,
    {
        let ((r, c), (_, cr)) = (self.0.shape_generic(), rhs.shape_generic());
        let mut out = OMatrix::zeros_generic(r, c.add(cr).sub(Const::<1>));
        self.mul_to(rhs, &mut out);
        Polynomial(out)
    }

    /// Mutliply `self` with `rhs` and store the result into `out`.
    pub fn mul_to<CR, SR, SO>(
        &self,
        rhs: &Matrix<T, R, CR, SR>,
        out: &mut Matrix<T, R, DimPolyProd<C, CR>, SO>,
    ) where
        T: Field,
        CR: Dim,
        C: DimAdd<CR>,
        DimSum<C, CR>: DimSub<U1>,
        SR: Storage<T, R, CR>,
        SO: StorageMut<T, R, DimPolyProd<C, CR>>,
    {
        out.fill(T::zero());
        for (i, x) in self.0.column_iter().enumerate() {
            for (j, y) in rhs.column_iter().enumerate() {
                for row in 0..self.0.nrows() {
                    out[(row, i + j)] += x[row].clone() * y[row].clone();
                }
            }
        }
    }
}
/// The column's dimension for the product of two polynomials
pub type DimPolyProd<C1, C2> = DimDiff<DimSum<C1, C2>, U1>;

/* Roots */
impl<T: Scalar, S: Storage<T, U1, U2>> Polynomial<T, U1, U2, S> {
    /// Calculate a linear's root.
    pub fn roots(&self) -> Vec<T>
    where
        T: RealField,
    {
        let t = &self.0[(0, 0)];
        let m = &self.0[(0, 1)];

        if m.is_zero() {
            vec![]
        } else {
            vec![(-t.clone()) / m.clone()]
        }
    }
}
impl<T: Scalar, S: Storage<T, U1, U3>> Polynomial<T, U1, U3, S> {
    /// Calculate a quadratic's roots.
    pub fn roots(&self) -> Vec<T>
    where
        T: RealField,
    {
        let two = T::one() + T::one();
        let four = two.clone() + two.clone();

        let c = &self.0[(0, 0)];
        let b = &self.0[(0, 1)];
        let a = &self.0[(0, 2)];

        let d: T = b.clone().powi(2) - four * a.clone() * c.clone();
        if d.is_sign_negative() {
            Vec::new()
        } else if d.is_zero() {
            vec![(-b.clone()) / (two * a.clone())]
        } else {
            let p = -b.clone();
            let q = d.sqrt();
            let r = two * a.clone();
            vec![(p.clone() + q.clone()) / r.clone(), (p - q) / r]
        }
    }
}
impl<T: Scalar, S: Storage<T, U1, U4>> Polynomial<T, U1, U4, S> {
    /// Calculate a cubic's roots.
    ///
    /// Based on [wikipedia](https://en.wikipedia.org/wiki/Cubic_equation#General_cubic_formula)'s formula.
    pub fn roots(&self) -> Vec<T>
    where
        T: RealField,
    {
        let two = T::one() + T::one();
        let three = T::one() + T::one() + T::one();
        let four = two.clone() + two.clone();
        let nine = three.clone() + three.clone() + three.clone();
        let twenty_seven = nine.clone() + nine.clone() + nine.clone();

        let d = &self.0[(0, 0)];
        let c = &self.0[(0, 1)];
        let b = &self.0[(0, 2)];
        let a = &self.0[(0, 3)];

        let delta_0 = b.clone().powi(2) - three.clone() * a.clone() * c.clone();
        let delta_1 = two.clone() * b.clone().powi(3)
            - nine.clone() * a.clone() * b.clone() * c.clone()
            + twenty_seven.clone() * a.clone().powi(2) * d.clone();

        if delta_0 == T::zero() && delta_1 == T::zero() {
            return vec![b.clone().neg() / three.clone() * a.clone()];
        }

        // Probably right up to here

        let pre_sqrt = delta_1.clone().powi(2) - four.clone() * delta_0.clone().powi(3);
        let sqrt = if pre_sqrt.is_sign_negative() {
            (T::zero(), pre_sqrt.neg().sqrt())
        } else {
            (pre_sqrt.sqrt(), T::zero())
        };

        let pre_cbrt = {
            let (real, imag) = sqrt;
            ((delta_1.clone() + real) / two.clone(), imag / two.clone())
        };
        let pre_cbrt = {
            let (real, imag) = pre_cbrt;
            (
                (real.clone().powi(2) + imag.clone().powi(2)).sqrt(),
                imag.atan2(real),
            )
        };
        let mut big_c = {
            let (r, phi) = pre_cbrt;
            (r.cbrt(), phi / three.clone())
        };

        let mut roots = Vec::new();
        for _ in 0..3 {
            let pre_x = (
                big_c.0.clone() * big_c.1.clone().cos()
                    + delta_0.clone() * big_c.0.clone().recip() * big_c.1.clone().neg().cos(),
                big_c.0.clone() * big_c.1.clone().sin()
                    + delta_0.clone() * big_c.0.clone().recip() * big_c.1.clone().neg().sin(),
            );

            if pre_x.1.abs() <= T::one().powi(-14) {
                let x = (three.clone() * a.clone()).recip().neg() * (b.clone() + pre_x.0);
                roots.push(x);
            }

            big_c.1 += T::two_pi() / three.clone();
        }

        roots
    }
}
impl<T: Scalar, S: Storage<T, U1, Dynamic>> Polynomial<T, U1, Dynamic, S> {
    /// Calculate `self`'s roots.
    ///
    /// **If `self` is any polynomial other than a quadratic or cubic one,
    /// this method will return an empty Vec, since only those two are supported yet.**
    pub fn roots(&self) -> Vec<T>
    where
        T: RealField,
    {
        match self.0.ncols() {
            2 => Polynomial(self.0.fixed_columns::<2>(0)).roots(),
            3 => Polynomial(self.0.fixed_columns::<3>(0)).roots(),
            4 => Polynomial(self.0.fixed_columns::<4>(0)).roots(),
            _ => vec![],
        }
    }
}

/* Common traits */
impl<T: Scalar, R: Dim, R2: Dim, C: Dim, C2: Dim, S, S2> PartialEq<Polynomial<T, R2, C2, S2>>
    for Polynomial<T, R, C, S>
where
    S: RawStorage<T, R, C>,
    S2: RawStorage<T, R2, C2>,
{
    #[inline]
    fn eq(&self, rhs: &Polynomial<T, R2, C2, S2>) -> bool {
        self.0 == rhs.0
    }
}
impl<T, R, C, S: fmt::Debug> fmt::Debug for Polynomial<T, R, C, S> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> Result<(), fmt::Error> {
        self.0.data.fmt(formatter)
    }
}

/* barely working Display implementation */
impl<T, R: Dim, C: Dim, S> fmt::Display for Polynomial<T, R, C, S>
where
    T: Scalar + fmt::Display,
    S: Storage<T, R, C>,
    DefaultAllocator: Allocator<T, R>,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let strings: Vec<_> = self
            .0
            .column_iter()
            .map(|c| format!("{}", c.into_owned()))
            .collect();

        let mut splits: Vec<_> = strings.iter().map(|s| s.split("\n")).collect();

        for row in 0..(self.0.nrows() + 5) {
            for col in 0..self.0.ncols() {
                if let Some(next) = splits[col].next() {
                    f.write_str(next)?;
                } else {
                    return Err(fmt::Error);
                }

                if row == 2 + self.0.nrows() / 2 {
                    write!(f, " x^{} +", col)?;
                } else {
                    let chars = 5 + if col < 10 { 1 } else { 2 };
                    for _ in 0..chars {
                        f.write_char(' ')?;
                    }
                }
            }
            f.write_char('\n')?;
        }
        Ok(())
    }
}
