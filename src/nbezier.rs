//! A wrapper around [`nalgebra::Matrix`] interpreting it as a bezier curve.

use nalgebra::allocator::Allocator;
use nalgebra::constraint::{DimEq, ShapeConstraint};
use nalgebra::dimension::{Const, Dim, DimAdd, DimDiff, DimSub, DimSum, U1, U2};
use nalgebra::{DefaultAllocator, Matrix, OMatrix, OVector, Owned, RealField, Storage, Vector2};

use crate::bounding_box::BoundingBox;
use crate::graham_scan::convex_hull;
use crate::npolynomial::Polynomial;

/// Wrapper around [`nalgebra::Matrix`] interpreting it as a bezier curve.
///
/// The curve's control points are stored as the matrix' columns.
pub struct BezierCurve<T, R, C, S>(pub Matrix<T, R, C, S>);

/// Wrapper around [`nalgebra::OMatrix`] interpreting it as a bezier curve.
pub type OBezierCurve<T, R, C> = BezierCurve<T, R, C, Owned<T, R, C>>;

impl<T: RealField, R: Dim, C: Dim, S: Storage<T, R, C>> BezierCurve<T, R, C, S> {
    /// Get the curves degree
    ///
    /// For example a cubic curve has degree 3 and 4 control points
    pub fn degree(&self) -> usize {
        self.0.ncols() - 1
    }
}

impl<T: RealField, R: Dim, C: Dim, S: Storage<T, R, C>> BezierCurve<T, R, C, S>
where
    // Column arithemtic required in each step
    DefaultAllocator: Allocator<T, R, U1>,

    // Buffer to store intermediate steps in
    DefaultAllocator: Allocator<T, R, C>,
{
    /// Splits a curve into two parts
    ///
    /// The first part is the same shape as the original curve between 0 and t and the second
    /// part as the curve between t and 1.
    /// This method assumes `t` to between 0 and 1 but doesn't check it.
    pub fn split(&self, t: T) -> [OBezierCurve<T, R, C>; 2] {
        let (rows, cols) = self.0.shape_generic();
        let mut lower: Matrix<T, R, C, _> = Matrix::zeros_generic(rows, cols);
        let mut upper: Matrix<T, R, C, _> = Matrix::zeros_generic(rows, cols);
        let ncols = cols.value();

        lower.set_column(0, &self.0.column(0));
        upper.set_column(ncols - 1, &self.0.column(ncols - 1));

        let t_inv = T::one() - t.clone();

        let mut points = (&mut self.0.clone_owned(), &mut self.0.clone_owned());
        for i in 1..ncols {
            let (input, output) = points;

            // castlejau step
            for (i, (p, q)) in input
                .column_iter()
                .zip(input.column_iter().skip(1))
                .enumerate()
            {
                let column = p * t_inv.clone() + q * t.clone();
                output.set_column(i, &column);
            }

            // copy edge points
            lower.set_column(i, &output.column(0));
            upper.set_column(ncols - 1 - i, &output.column(ncols - i - 1));

            points = (output, input);
        }

        [BezierCurve(lower), BezierCurve(upper)]
    }

    /// Get the point on the curve at position `t`.
    ///
    /// This method uses de castlejau's algorithm. An alternative way would be to evaluate the
    /// curve's polynomial (See `BezierCurve::polynomial`).
    pub fn castlejau_eval(&self, t: T) -> OVector<T, R> {
        let t_inv = T::one() - t.clone();
        let ncols = self.0.ncols();

        let mut points = (&mut self.0.clone_owned(), &mut self.0.clone_owned());
        for step in 0..ncols {
            let (input, output) = points;
            for i in 1..(ncols - step) {
                let column = &input.column(i - 1) * t_inv.clone() + &input.column(i) * t.clone();
                output.set_column(i - 1, &column);
            }
            points = (output, input);
        }

        return points.1.column(0).clone_owned();
    }
}

impl<T: RealField, R: Dim, C: Dim, S: Storage<T, R, C>> BezierCurve<T, R, C, S>
where
    C: DimAdd<U1>,
    DefaultAllocator: Allocator<T, C, R>,             // self^T
    DefaultAllocator: Allocator<T, DimSum<C, U1>, C>, // M
    DefaultAllocator: Allocator<T, DimSum<C, U1>, R>, // result^T
    DefaultAllocator: Allocator<T, R, DimSum<C, U1>>, // result
{
    /// Returns a new bezier curve of the same shape whose degree is one step higher
    pub fn raise_degree(&self) -> OBezierCurve<T, R, DimSum<C, U1>> {
        let (_, cols) = self.0.shape_generic();
        let matrix = elevation_matrix(cols);

        // Apply transformation
        let old = self.0.transpose();
        let new = matrix * &old;
        BezierCurve(new.transpose())
    }
}

impl<T: RealField, R: Dim, C: Dim, S: Storage<T, R, C>> BezierCurve<T, R, C, S>
where
    C: DimSub<U1>,
    DimDiff<C, U1>: DimAdd<U1>,
    ShapeConstraint: DimEq<DimSum<DimDiff<C, U1>, U1>, C>,
    DefaultAllocator: Allocator<T, DimSum<DimDiff<C, U1>, U1>, DimDiff<C, U1>>, // M
    DefaultAllocator: Allocator<T, DimDiff<C, U1>, DimSum<DimDiff<C, U1>, U1>>, // M^T
    DefaultAllocator: Allocator<T, DimDiff<C, U1>, DimDiff<C, U1>>,             // M^T M
    DefaultAllocator: Allocator<T, C, R>,                                       // self^T
    DefaultAllocator: Allocator<T, DimDiff<C, U1>, R>,                          // result^T
    DefaultAllocator: Allocator<T, R, DimDiff<C, U1>>,                          // result
{
    /// Returns a new bezier curve of an approximated shape whose degree is one step lower
    pub fn reduce_degree(&self) -> OBezierCurve<T, R, DimDiff<C, U1>> {
        let (_, cols) = self.0.shape_generic();

        // Compute (M^T M)^{-1} M^T
        let matrix = elevation_matrix(cols.sub(Const::<1>));
        let transposed = matrix.transpose();
        let square = matrix.tr_mul(&matrix);
        let inverse = square.try_inverse().unwrap();
        let matrix = inverse * transposed;

        // Apply transformation
        let old = self.0.transpose();
        let new = matrix * old;
        BezierCurve(new.transpose())
    }
}

/// Constructs the matrix to raise a curve's degree from `n` to `n+1`
fn elevation_matrix<T: RealField, N: Dim>(n: N) -> OMatrix<T, DimSum<N, U1>, N>
where
    N: DimAdd<U1>,
    DefaultAllocator: Allocator<T, DimSum<N, U1>, N>,
{
    let mut matrix = OMatrix::zeros_generic(n.add(Const::<1>), n);

    let n = n.value();
    let n_as_t = {
        let mut k = T::zero();
        for _ in 0..n {
            k = k + T::one();
        }
        k
    };

    matrix[(0, 0)] = T::one();
    matrix[(n, n - 1)] = T::one();

    let mut i_as_t = T::zero();
    for i in 1..n {
        i_as_t += T::one();
        matrix[(i, i - 1)] = i_as_t.clone() / n_as_t.clone();
        matrix[(i, i)] = T::one() - matrix[(i, i - 1)].clone();
    }
    matrix
}

impl<T: RealField, R: Dim, C: Dim, S: Storage<T, R, C>> BezierCurve<T, R, C, S>
where
    DefaultAllocator: Allocator<T, R, C>, // polynomial
    DefaultAllocator: Allocator<T, C, C>, // bernstein basis
{
    /// Computes the curve's polynomial
    ///
    /// This polynomial evaluated between `0` and `1` yields the same points as its corrisponding bezier curve.
    ///
    /// If you are only interested in its derivative, use [`derivative`] to get it directly.
    ///
    /// [`derivative`]: BezierCurve::derivative
    pub fn polynomial(&self) -> Polynomial<T, R, C, Owned<T, R, C>> {
        let (rows, cols) = self.0.shape_generic();
        let mut polynomial = Matrix::zeros_generic(rows, cols);

        for (i, b) in bernstein_polynomials::<T, _>(cols).row_iter().enumerate() {
            let p = self.0.column(i) * b;
            polynomial += p;
        }

        Polynomial(polynomial)
    }
}

impl<T: RealField, R: Dim, C: Dim, S: Storage<T, R, C>> BezierCurve<T, R, C, S>
where
    C: DimSub<U1>,
    DefaultAllocator: Allocator<T, R, DimDiff<C, U1>>, // polynomial
    DefaultAllocator: Allocator<T, DimDiff<C, U1>, DimDiff<C, U1>>, // bernstein basis
    DefaultAllocator: Allocator<T, R, U1>, // column difference TODO: might be avoidable?
{
    /// Computes the curve's polynomial's derivative
    ///
    /// This method is a faster alternative to calling [`Polynomial::derive`] on the result of
    /// [`BezierCurve::polynomial`].
    pub fn derivative(&self) -> Polynomial<T, R, DimDiff<C, U1>, Owned<T, R, DimDiff<C, U1>>> {
        let (rows, cols) = self.0.shape_generic();
        let cols = cols.sub(Const::<1>);
        let mut polynomial = Matrix::zeros_generic(rows, cols);

        for (i, b) in bernstein_polynomials::<T, _>(cols).row_iter().enumerate() {
            let p = (self.0.column(i + 1) - self.0.column(i)) * b;
            polynomial += p;
        }

        let degree = cols.value();
        let degree = {
            let mut k = T::zero();
            for _ in 0..degree {
                k += T::one();
            }
            k
        };
        polynomial *= degree;

        Polynomial(polynomial)
    }

    /// Computes the curve's tangent vector at `t`
    ///
    /// If you need a lot of them at once, it is far more efficent to call
    /// [`derivative`] once yourself and evaluate it multiple times,
    /// as this function would recompute the derivative for every vector.
    ///
    /// *The resulting vector is not normalized!*
    ///
    /// [`derivative`]: BezierCurve::derivative
    pub fn tangent(&self, t: T) -> OVector<T, R> {
        self.derivative().evaluate(t)
    }
}

impl<T: RealField, C: Dim, S: Storage<T, U2, C>> BezierCurve<T, U2, C, S>
where
    C: DimSub<U1>,
    DefaultAllocator: Allocator<T, U2, DimDiff<C, U1>>, // polynomial
    DefaultAllocator: Allocator<T, DimDiff<C, U1>, DimDiff<C, U1>>, // bernstein basis
    DefaultAllocator: Allocator<T, U2, U1>, // column difference TODO: might be avoidable?
{
    /// Computes the curve's normal vector at `t`
    ///
    /// Similarly to [`tangent`] calling this multiple times to get a lot of vectors
    /// would be inefficent. Compute their tangent vectors (see [`tangent`]) and rotate
    /// them yourself instead.
    ///
    /// *The resulting vector is not normalized!*
    ///
    /// [`tangent`]: BezierCurve::tangent
    pub fn normal(&self, t: T) -> Vector2<T> {
        let tangent = self.tangent(t);
        Vector2::new(T::zero() - tangent.y.clone(), tangent.x.clone())
    }
}

/// Computes the bernstein polynomial basis for a given degree
pub fn bernstein_polynomials<T: RealField, C: Dim>(cols: C) -> OMatrix<T, C, C>
where
    DefaultAllocator: Allocator<T, C, C>,
{
    // Each row is a different berstein polynomial
    let mut polynomials = OMatrix::zeros_generic(cols, cols);

    // Fill matrix with pascal triangle of shape:
    //    ...
    //   1 3 3 1
    //     1 2 1
    //       1 1
    //         1
    for n in 0..cols.value() {
        let m = cols.value() - 1 - n;
        polynomials[(m, cols.value() - 1)] = T::one();
        polynomials[(m, m)] = T::one();
        for i in 1..n {
            let j = cols.value() - 1 - i;
            polynomials[(m, j)] =
                polynomials[(m + 1, j + 1)].clone() + polynomials[(m + 1, j)].clone();
        }
    }

    // Scale every row by the entries is the top row
    // (First and last row can be skipped, since their coeff is always 1)
    for i in 1..(cols.value() - 1) {
        let coeff = polynomials[(0, i)].clone();
        let mut row = polynomials.row_mut(i);
        row *= coeff;
    }

    // Apply minus sign in checkerboard pattern
    for i in 0..cols.value() {
        for j in i..cols.value() {
            if (j + i) % 2 == 1 {
                polynomials[(i, j)] *= T::one().neg();
            }
        }
    }

    polynomials
}

impl<T: RealField, C: Dim, S: Storage<T, U2, C>> BezierCurve<T, U2, C, S> {
    /// Constructs an axis aligned bounding box containing all controll points.
    ///
    /// This box will also contain the whole curve, but can highly overestimate it.
    /// It can be used as a fast way to estimate intersections.
    /// For a more precise boundary consider: [`convex_hull`]
    ///
    /// [`convex_hull`]: BezierCurve::convex_hull
    pub fn bounding_box(&self) -> BoundingBox<T> {
        BoundingBox::from_iter(self.0.column_iter().map(|column| column.clone_owned()))
    }

    /// Computes the control points' convex hull using [graham scan](https://en.wikipedia.org/wiki/Graham_scan).
    ///
    /// *Actual code for intersecting convex polygons is not implemented yet!*
    pub fn convex_hull(&self) -> Vec<Vector2<T>> {
        convex_hull(
            self.0
                .column_iter()
                .map(|column| column.clone_owned())
                .collect(),
        )
    }

    //TODO Minimal Bounding Box isn't possible due to non generic constraints to Polynomial::roots

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
    pub fn locate_point(&self, p: Vector2<T>) -> Option<T>
    where
        // Clone curves
        DefaultAllocator: Allocator<T, U2, C>,

        // Tangent
        C: DimSub<U1>,
        DefaultAllocator: Allocator<T, U2, DimDiff<C, U1>>,
        DefaultAllocator: Allocator<T, DimDiff<C, U1>, DimDiff<C, U1>>,
    {
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
        let mut divisions = (
            &mut vec![SubCurve {
                from: T::zero(),
                to: T::one(),
                curve: BezierCurve(self.0.clone_owned()),
            }],
            &mut Vec::new(),
        );
        for _ in 0..NUM_SUBDIVISIONS {
            if divisions.0.len() == 0 {
                return None;
            }
            for subcurve in divisions.0.iter() {
                let bb = subcurve.curve.bounding_box();
                if bb.contains(p.clone()) {
                    let [lower, upper] = subcurve.split();
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
    pub fn get_intersections<CO, SO>(&self, other: &BezierCurve<T, U2, CO, SO>) -> Vec<Vector2<T>>
    where
        CO: Dim,
        SO: Storage<T, U2, CO>,
        // Clone curves
        DefaultAllocator: Allocator<T, U2, C>,
        DefaultAllocator: Allocator<T, U2, CO>,
    {
        #[allow(non_snake_case)]
        let NUM_SUBDIVISIONS: usize = 20;

        // Subdivide curves and check bounding boxes
        let mut divisions = (
            &mut vec![(
                SubCurve {
                    from: T::zero(),
                    to: T::one(),
                    curve: BezierCurve(self.0.clone_owned()),
                },
                SubCurve {
                    from: T::zero(),
                    to: T::one(),
                    curve: BezierCurve(other.0.clone_owned()),
                },
            )],
            &mut Vec::new(),
        );
        for _ in 0..NUM_SUBDIVISIONS {
            if divisions.0.len() == 0 {
                return Vec::new();
            }
            for (s, o) in divisions.0.iter() {
                if s.curve.bounding_box().intersects(&o.curve.bounding_box()) {
                    let [s_lower, s_upper] = s.split();
                    let [o_lower, o_upper] = o.split();
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

/* Helper struct repeatedly subdividing a curve */
struct SubCurve<T, R, C, S> {
    from: T,
    to: T,
    curve: BezierCurve<T, R, C, S>,
}
impl<T: RealField, C: Dim, S: Storage<T, U2, C>> SubCurve<T, U2, C, S>
where
    DefaultAllocator: Allocator<T, U2, C>,
    DefaultAllocator: Allocator<T, U2, U1>,
{
    fn split(&self) -> [SubCurve<T, U2, C, Owned<T, U2, C>>; 2] {
        let two = T::one() + T::one();
        let middle = (self.from.clone() + self.to.clone()) / two.clone();

        let [lower, upper] = self.curve.split(T::one() / two);
        [
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
        ]
    }

    fn middle_point(&self) -> OVector<T, U2> {
        let two = T::one() + T::one();
        let middle = (self.from.clone() + self.to.clone()) / two;

        self.curve.castlejau_eval(middle)
    }
}
impl<T: RealField, R: Dim, C: Dim> Clone for SubCurve<T, R, C, Owned<T, R, C>>
where
    DefaultAllocator: Allocator<T, R, C>,
{
    fn clone(&self) -> Self {
        SubCurve {
            from: self.from.clone(),
            to: self.to.clone(),
            curve: BezierCurve(self.curve.0.clone_owned()),
        }
    }
}
