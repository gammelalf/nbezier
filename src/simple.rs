//! Wrapper for [BezierCurve] and [Polynomial] with simple interface
//!
//! These wrapper reintroduce this library's original idea:
//! A single datatype storing curves of arbirtary degree which is optimised for cubic curves and lower.

use crate::bounding_box::BoundingBox;
use crate::nbezier::OBezierCurve;
use crate::npolynomial::{OPolynomial, Polynomial};
use crate::BezierCurve;
use nalgebra::allocator::Allocator;
use nalgebra::default_allocator::DefaultAllocator;
use nalgebra::dimension::{Const, Dim, Dynamic};
use nalgebra::{OMatrix, Vector2};

/// Bezier Curve of arbitrary degree which stores cubic curves and lower degrees on the stack.
pub enum SimpleCurve {
    /// Stack allocated linear curve i.e. a straight line
    Linear(OBezierCurve<f64, Const<2>, Const<2>>),

    /// Stack allocated quadratic curve
    Quadratic(OBezierCurve<f64, Const<2>, Const<3>>),

    /// Stack allocated cubic curve
    Cubic(OBezierCurve<f64, Const<2>, Const<4>>),

    /// Heap allocated curve of arbitrary degree
    Higher(OBezierCurve<f64, Const<2>, Dynamic>),
}

impl SimpleCurve {
    /// Splits a curve into two parts
    ///
    /// The first part is the same shape as the original curve between 0 and t and the second
    /// part as the curve between t and 1.
    /// This method assumes `t` to between 0 and 1 but doesn't check it.
    pub fn split(&self, t: f64) -> [Self; 2] {
        match self {
            SimpleCurve::Linear(curve) => curve.split(t).map(From::from),
            SimpleCurve::Quadratic(curve) => curve.split(t).map(From::from),
            SimpleCurve::Cubic(curve) => curve.split(t).map(From::from),
            SimpleCurve::Higher(curve) => curve.split(t).map(From::from),
        }
    }

    /// Get the point on the curve at position `t`.
    ///
    /// This method uses de castlejau's algorithm. An alternative way would be to evaluate the
    /// curve's polynomial (See `SimpleCurve::polynomial`).
    pub fn castlejau_eval(&self, t: f64) -> Vector2<f64> {
        match self {
            SimpleCurve::Linear(curve) => curve.castlejau_eval(t),
            SimpleCurve::Quadratic(curve) => curve.castlejau_eval(t),
            SimpleCurve::Cubic(curve) => curve.castlejau_eval(t),
            SimpleCurve::Higher(curve) => curve.castlejau_eval(t),
        }
    }

    /// Returns a new bezier curve of the same shape whose degree is one step higher
    pub fn raise_degree(&self) -> Self {
        match self {
            SimpleCurve::Linear(curve) => curve.raise_degree().into(),
            SimpleCurve::Quadratic(curve) => curve.raise_degree().into(),
            SimpleCurve::Cubic(curve) => {
                let matrix = curve.raise_degree().0.resize_horizontally(0, 0.0);
                SimpleCurve::Higher(BezierCurve(matrix))
            }
            SimpleCurve::Higher(curve) => curve.raise_degree().into(),
        }
    }

    /// Returns a new bezier curve of an approximated shape whose degree is one step lower
    pub fn reduce_degree(&self) -> Self {
        match self {
            SimpleCurve::Linear(curve) => SimpleCurve::Linear(BezierCurve(curve.0.clone_owned())),
            SimpleCurve::Quadratic(curve) => curve.reduce_degree().into(),
            SimpleCurve::Cubic(curve) => curve.reduce_degree().into(),
            SimpleCurve::Higher(curve) => curve.reduce_degree().into(),
        }
    }

    /// Computes the curve's polynomial
    ///
    /// This polynomial evaluated between `0` and `1` yields the same points as its corrisponding bezier curve.
    ///
    /// If you are only interested in its derivative, use [`derivative`] to get it directly.
    ///
    /// [`derivative`]: SimpleCurve::derivative
    pub fn polynomial(&self) -> SimplePolynomial {
        match self {
            SimpleCurve::Linear(curve) => curve.polynomial().into(),
            SimpleCurve::Quadratic(curve) => curve.polynomial().into(),
            SimpleCurve::Cubic(curve) => curve.polynomial().into(),
            SimpleCurve::Higher(curve) => curve.polynomial().into(),
        }
    }

    /// Computes the curve's polynomial's derivative
    ///
    /// This method is a faster alternative to calling [`SimplePolynomial::derive`] on the result of
    /// [`SimpleCurve::polynomial`].
    pub fn derivative(&self) -> SimplePolynomial {
        match self {
            SimpleCurve::Linear(curve) => curve.derivative().into(),
            SimpleCurve::Quadratic(curve) => curve.derivative().into(),
            SimpleCurve::Cubic(curve) => curve.derivative().into(),
            SimpleCurve::Higher(curve) => curve.derivative().into(),
        }
    }

    /// Computes the curve's tangent vector at `t`
    ///
    /// If you need a lot of them at once, it is far more efficent to call
    /// [`derivative`] once yourself and evaluate it multiple times,
    /// as this function would recompute the derivative for every vector.
    ///
    /// *The resulting vector is not normalized!*
    ///
    /// [`derivative`]: SimpleCurve::derivative
    pub fn tangent(&self, t: f64) -> Vector2<f64> {
        self.derivative().evaluate(t)
    }

    /// Computes the curve's normal vector at `t`
    ///
    /// Similarly to [`tangent`] calling this multiple times to get a lot of vectors
    /// would be inefficent. Compute their tangent vectors (see [`tangent`]) and rotate
    /// them yourself instead.
    ///
    /// *The resulting vector is not normalized!*
    ///
    /// [`tangent`]: SimpleCurve::tangent
    pub fn normal(&self, t: f64) -> Vector2<f64> {
        let tangent = self.tangent(t);
        Vector2::new(0.0 - tangent.y, tangent.x)
    }

    /// Constructs an axis aligned bounding box containing all controll points.
    ///
    /// This box will also contain the whole curve, but can highly overestimate it.
    /// It can be used as a fast way to estimate intersections.
    /// For a more precise boundary consider: [`convex_hull`]
    ///
    /// [`convex_hull`]: SimpleCurve::convex_hull
    pub fn bounding_box(&self) -> BoundingBox<f64> {
        match self {
            SimpleCurve::Linear(curve) => curve.bounding_box(),
            SimpleCurve::Quadratic(curve) => curve.bounding_box(),
            SimpleCurve::Cubic(curve) => curve.bounding_box(),
            SimpleCurve::Higher(curve) => curve.bounding_box(),
        }
    }

    /// Computes the control points' convex hull using [graham scan](https://en.wikipedia.org/wiki/Graham_scan).
    ///
    /// *Actual code for intersecting convex polygons is not implemented yet!*
    pub fn convex_hull(&self) -> Vec<Vector2<f64>> {
        match self {
            SimpleCurve::Linear(curve) => curve.convex_hull(),
            SimpleCurve::Quadratic(curve) => curve.convex_hull(),
            SimpleCurve::Cubic(curve) => curve.convex_hull(),
            SimpleCurve::Higher(curve) => curve.convex_hull(),
        }
    }

    /// Trys to locate a point on the curve and returns its `t` value if located.
    pub fn locate_point(&self, p: Vector2<f64>) -> Option<f64> {
        match self {
            SimpleCurve::Linear(curve) => curve.locate_point(p),
            SimpleCurve::Quadratic(curve) => curve.locate_point(p),
            SimpleCurve::Cubic(curve) => curve.locate_point(p),
            SimpleCurve::Higher(curve) => curve.locate_point(p),
        }
    }
}

impl From<OBezierCurve<f64, Const<2>, Const<2>>> for SimpleCurve {
    /// Wrap a linear curve into a simple one
    fn from(curve: OBezierCurve<f64, Const<2>, Const<2>>) -> Self {
        SimpleCurve::Linear(curve)
    }
}
impl From<OBezierCurve<f64, Const<2>, Const<3>>> for SimpleCurve {
    /// Wrap a quadratic curve into a simple one
    fn from(curve: OBezierCurve<f64, Const<2>, Const<3>>) -> Self {
        SimpleCurve::Quadratic(curve)
    }
}
impl From<OBezierCurve<f64, Const<2>, Const<4>>> for SimpleCurve {
    /// Wrap a cubic curve into a simple one
    fn from(curve: OBezierCurve<f64, Const<2>, Const<4>>) -> Self {
        SimpleCurve::Cubic(curve)
    }
}
impl From<OBezierCurve<f64, Const<2>, Dynamic>> for SimpleCurve {
    /// Wrap a dynamic curve into a simple one
    ///
    /// If the dynamic curve is actualy of degree 3 or lower,
    /// it will be converted into a curve of the apropriate static degree.
    fn from(curve: OBezierCurve<f64, Const<2>, Dynamic>) -> Self {
        match curve.degree() {
            1 => {
                let matrix = curve.0.columns_generic(0, Const::<2>).clone_owned();
                SimpleCurve::Linear(BezierCurve(matrix))
            }
            2 => {
                let matrix = curve.0.columns_generic(0, Const::<3>).clone_owned();
                SimpleCurve::Quadratic(BezierCurve(matrix))
            }
            3 => {
                let matrix = curve.0.columns_generic(0, Const::<4>).clone_owned();
                SimpleCurve::Cubic(BezierCurve(matrix))
            }
            _ => SimpleCurve::Higher(curve),
        }
    }
}
impl<C: Dim, R: Dim> From<OMatrix<f64, C, R>> for SimpleCurve
where
    SimpleCurve: From<OBezierCurve<f64, C, R>>,
    DefaultAllocator: Allocator<f64, C, R>,
{
    /// Wrap any matrix, whose curve is wrappable, directly into a simple curve
    fn from(matrix: OMatrix<f64, C, R>) -> Self {
        BezierCurve(matrix).into()
    }
}

/// Polynomial of arbitrary degree which stores cubic polynomials and lower degrees on the stack.
pub enum SimplePolynomial {
    /// Stack allocated constant polynomial
    Constant(OPolynomial<f64, Const<2>, Const<1>>),

    /// Stack allocated linear polynomial i.e. `ax + b`
    Linear(OPolynomial<f64, Const<2>, Const<2>>),

    /// Stack allocated quadratic polynomial i.e. `ax^2 + bx + c`
    Quadratic(OPolynomial<f64, Const<2>, Const<3>>),

    /// Stack allocated cubic polynomial i.e. `ax^3 + bx^2 + cx + d`
    Cubic(OPolynomial<f64, Const<2>, Const<4>>),

    /// Heap allocated curve of arbitrary degree
    Higher(OPolynomial<f64, Const<2>, Dynamic>),
}

impl SimplePolynomial {
    /// Evaluate `self` at position `x`.
    pub fn evaluate(&self, x: f64) -> Vector2<f64> {
        match self {
            SimplePolynomial::Constant(poly) => poly.0.clone_owned(),
            SimplePolynomial::Linear(poly) => poly.evaluate(x),
            SimplePolynomial::Quadratic(poly) => poly.evaluate(x),
            SimplePolynomial::Cubic(poly) => poly.evaluate(x),
            SimplePolynomial::Higher(poly) => poly.evaluate(x),
        }
    }

    /// Calculate `self`'s derivative.
    pub fn derive(&self) -> SimplePolynomial {
        match self {
            SimplePolynomial::Constant(_) => Polynomial(Vector2::zeros()).into(),
            SimplePolynomial::Linear(poly) => poly.derive().into(),
            SimplePolynomial::Quadratic(poly) => poly.derive().into(),
            SimplePolynomial::Cubic(poly) => poly.derive().into(),
            SimplePolynomial::Higher(poly) => poly.derive().into(),
        }
    }

    /// Calculate `self`'s integral.
    ///
    /// The zero vector is used as integration constant.
    pub fn integrate(&self) -> SimplePolynomial {
        match self {
            SimplePolynomial::Constant(poly) => poly.integrate().into(),
            SimplePolynomial::Linear(poly) => poly.integrate().into(),
            SimplePolynomial::Quadratic(poly) => poly.integrate().into(),
            SimplePolynomial::Cubic(poly) => poly.integrate().0.resize_horizontally(0, 0.0).into(),
            SimplePolynomial::Higher(poly) => poly.integrate().into(),
        }
    }
}

impl From<OPolynomial<f64, Const<2>, Const<1>>> for SimplePolynomial {
    /// Wrap a constant polynomial into a simple one
    fn from(polynomial: OPolynomial<f64, Const<2>, Const<1>>) -> Self {
        SimplePolynomial::Constant(polynomial)
    }
}
impl From<OPolynomial<f64, Const<2>, Const<2>>> for SimplePolynomial {
    /// Wrap a linear polynomial into a simple one
    fn from(polynomial: OPolynomial<f64, Const<2>, Const<2>>) -> Self {
        SimplePolynomial::Linear(polynomial)
    }
}
impl From<OPolynomial<f64, Const<2>, Const<3>>> for SimplePolynomial {
    /// Wrap a quadratic polynomial into a simple one
    fn from(polynomial: OPolynomial<f64, Const<2>, Const<3>>) -> Self {
        SimplePolynomial::Quadratic(polynomial)
    }
}
impl From<OPolynomial<f64, Const<2>, Const<4>>> for SimplePolynomial {
    /// Wrap a cubic polynomial into a simple one
    fn from(polynomial: OPolynomial<f64, Const<2>, Const<4>>) -> Self {
        SimplePolynomial::Cubic(polynomial)
    }
}
impl From<OPolynomial<f64, Const<2>, Dynamic>> for SimplePolynomial {
    /// Wrap a dynamic polynomial into a simple one
    ///
    /// If the dynamic polynomial is actualy of degree 3 or lower,
    /// it will be converted into a polynomial of the apropriate static degree.
    fn from(polynomial: OPolynomial<f64, Const<2>, Dynamic>) -> Self {
        match polynomial.degree() {
            1 => {
                let matrix = polynomial.0.fixed_columns::<2>(0).clone_owned();
                SimplePolynomial::Linear(Polynomial(matrix))
            }
            2 => {
                let matrix = polynomial.0.fixed_columns::<3>(0).clone_owned();
                SimplePolynomial::Quadratic(Polynomial(matrix))
            }
            3 => {
                let matrix = polynomial.0.fixed_columns::<4>(0).clone_owned();
                SimplePolynomial::Cubic(Polynomial(matrix))
            }
            _ => SimplePolynomial::Higher(polynomial),
        }
    }
}
impl<C: Dim, R: Dim> From<OMatrix<f64, C, R>> for SimplePolynomial
where
    SimplePolynomial: From<OPolynomial<f64, C, R>>,
    DefaultAllocator: Allocator<f64, C, R>,
{
    /// Wrap any matrix, whose polynomial is wrappable, directly into a simple polynomial
    fn from(matrix: OMatrix<f64, C, R>) -> Self {
        Polynomial(matrix).into()
    }
}
