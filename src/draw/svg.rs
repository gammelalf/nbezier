//! Draw a curve in a svg

use crate::draw::DrawCurve;
use crate::{BezierCurve, SimpleCurve};
use std::fmt::Write;
use std::marker::PhantomData;

/// Helper trait implemented by [`Absolute`] and [`Relative`] to choose
/// how the coordinates are interpreted in the svg's path.
///
/// This effectively toggles between lower and upper case letters for the path commands.
pub trait CoordRepr {
    /// Character to use for the MoveTo command
    const M: &'static str;

    /// Character to use for the LineTo command
    const L: &'static str;

    /// Character to use for the Quadratic Bezier Curve command
    const Q: &'static str;

    /// Character to use for the Cubic Bezier Curve command
    const C: &'static str;
}

/// Interpret the svg path's coordinates as absolute.
pub struct Absolute;
impl CoordRepr for Absolute {
    const M: &'static str = "M";
    const L: &'static str = "L";
    const Q: &'static str = "Q";
    const C: &'static str = "C";
}

/// Interpret the svg path's coordinates as relative.
pub struct Relative;
impl CoordRepr for Relative {
    const M: &'static str = "m";
    const L: &'static str = "l";
    const Q: &'static str = "q";
    const C: &'static str = "c";
}

/// This type wraps a mutable String reference and implements [`DrawCurve`] on it.
/// The [`DrawCurve::add_curve`] writes the curve to the string
/// using the path commands used in a svg `<path>`'s d attribute
///
/// ```
/// # use nalgebra::Matrix2;
/// # use nbezier::SimpleCurve;
/// use nbezier::draw::DrawCurve;
/// use nbezier::draw::svg::SVGAbsolutePath;
///
/// let curve = SimpleCurve::from(Matrix2::new(0.0, 1.0, 2.0, 3.0));
/// let mut d = String::new();
/// SVGAbsolutePath::from(&mut d).add_curve(&curve);
/// assert!_eq(d, "M 0,1 L 2,3");
/// ```
pub struct SVGPath<'s, R: CoordRepr>(&'s mut String, PhantomData<&'s R>);

/// Wrapper for writing a curve to a svg's path using absolute coordinates
pub type SVGAbsolutePath<'s> = SVGPath<'s, Absolute>;

/// Wrapper for writing a curve to a svg's path using relative coordinates
pub type SVGRelativePath<'s> = SVGPath<'s, Relative>;

impl<'s, R: CoordRepr> From<&'s mut String> for SVGPath<'s, R> {
    fn from(string: &'s mut String) -> Self {
        SVGPath(string, PhantomData)
    }
}

impl<'s, R: CoordRepr> DrawCurve for SVGPath<'s, R> {
    fn add_curve(&mut self, curve: &SimpleCurve) {
        match curve {
            SimpleCurve::Linear(BezierCurve(matrix)) => {
                let _ = write!(
                    self.0,
                    "{} {},{} {} {},{}",
                    R::M,
                    matrix.column(0).x,
                    matrix.column(0).y,
                    R::L,
                    matrix.column(1).x,
                    matrix.column(1).y
                );
            }
            SimpleCurve::Quadratic(BezierCurve(matrix)) => {
                let _ = write!(
                    self.0,
                    "{} {},{} {} {},{} {},{}",
                    R::M,
                    matrix.column(0).x,
                    matrix.column(0).y,
                    R::Q,
                    matrix.column(1).x,
                    matrix.column(1).y,
                    matrix.column(2).x,
                    matrix.column(2).y,
                );
            }
            SimpleCurve::Cubic(BezierCurve(matrix)) => {
                let _ = write!(
                    self.0,
                    "{} {},{} {} {},{} {},{} {},{}",
                    R::M,
                    matrix.column(0).x,
                    matrix.column(0).y,
                    R::C,
                    matrix.column(1).x,
                    matrix.column(1).y,
                    matrix.column(2).x,
                    matrix.column(2).y,
                    matrix.column(3).x,
                    matrix.column(3).y,
                );
            }
            SimpleCurve::Higher(curve) => {
                const STEPS: usize = 20;

                let matrix = &curve.0;
                let n = matrix.ncols();

                let _ = write!(
                    self.0,
                    "{} {},{}",
                    R::M,
                    matrix.column(0).x,
                    matrix.column(0).y,
                );
                for i in 1..STEPS {
                    let t = i as f64 / STEPS as f64;
                    let point = curve.castlejau_eval(t);
                    let _ = write!(self.0, "{} {},{}", R::L, point.x, point.y);
                }
                let _ = write!(
                    self.0,
                    "{} {},{}",
                    R::L,
                    matrix.column(n - 1).x,
                    matrix.column(n - 1).y,
                );
            }
        }
    }
}
