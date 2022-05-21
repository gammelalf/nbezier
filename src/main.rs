use crate::bezier::BezierCurve;
use crate::svg::SVG;
use crate::vector::Vector;

mod vector;
mod graham_scan;
mod polynomial;
mod bezier;
mod svg;

fn main() {
    let mut svg = SVG {
        view_box: (0.0, 0.0, 100.0, 100.0),
        elements: Vec::with_capacity(0),
    };
    let curve = BezierCurve(vec![
        Vector([50.0, 0.0]),
        Vector([200.0, 33.0]),
        Vector([0.0, 66.0]),
        Vector([50.0, 100.0]),
    ]);
    let (upper, lower) = curve.split(0.7).unwrap();
    svg.debug_bezier(&upper, "blue");
    svg.debug_bezier(&lower, "red");
    //svg.debug_bezier(&curve, "green");
    println!("{}", svg);
}