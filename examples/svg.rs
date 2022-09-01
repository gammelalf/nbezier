use nalgebra::{Matrix2xX, Vector2};
use nbezier::nbezier::BezierCurve;

#[path = "../src/svg.rs"]
mod svg;
use svg::SVG;

fn main() {
    let mut svg = SVG {
        view_box: (0.0, 0.0, 100.0, 100.0),
        elements: Vec::with_capacity(0),
    };
    let curve = BezierCurve(Matrix2xX::from_columns(&[
        Vector2::new(50.0, 0.0),
        Vector2::new(200.0, 33.0),
        Vector2::new(0.0, 66.0),
        Vector2::new(50.0, 100.0),
    ]));

    let [lower, upper] = curve.split(0.7);
    svg.debug_bezier(&upper, "red");
    svg.debug_bezier(&lower, "blue");

    println!("{}", svg);
}
