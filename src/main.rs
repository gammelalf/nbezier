use crate::bezier::BezierCurve;
use crate::svg::SVG;
use nalgebra::Vector2;
use smallvec::smallvec;

mod bezier;
mod bounding_box;
mod graham_scan;
mod npolynomial;
mod svg;

fn main() {
    let mut svg = SVG {
        view_box: (0.0, 0.0, 100.0, 100.0),
        elements: Vec::with_capacity(0),
    };
    let curve = BezierCurve(smallvec![
        Vector2::new(50.0, 0.0),
        Vector2::new(200.0, 33.0),
        Vector2::new(0.0, 66.0),
        Vector2::new(50.0, 100.0),
    ]);
    let (upper, lower) = curve.split(0.7);
    svg.debug_bezier(&upper, "blue");
    svg.debug_bezier(&lower, "red");
    println!("{}", svg);
}
