//! Defines and implements the helper method [`DrawCurve::add_curve`] on various drawing contexts.
//!
//! Drawing a curve (especially a cubic one which most drawing contexts support directly)
//! is never really an issue.
//! But it is noisy, repetitive and annoyed me across various test projects.
use crate::SimpleCurve;

/// Draw bezier curves on different "drawing contexts" with ease.
///
/// Use different crate features to implement different contexts:
/// - `draw-svg` to draw curves using svg paths
pub trait DrawCurve {
    /// Add a curve to the drawing context.
    ///
    /// What this actually means depends on the context.
    /// But generally this just prepares the curve to be drawn instead of actually drawing it.
    fn add_curve(&mut self, curve: &SimpleCurve);
}

#[cfg(feature = "draw-svg")]
pub mod svg;
