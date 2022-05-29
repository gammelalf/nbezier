use serde::Deserialize;
use wasm_bindgen::prelude::*;
use js_sys::Array;
use web_sys::CanvasRenderingContext2d;
use gammalg::vector::Vector;
use gammalg::bezier::BezierCurve;
use gammalg::graham_scan::convex_hull;

#[derive(Deserialize)]
pub struct Point {
    pub x: f64,
    pub y: f64,
}
impl From<Point> for Vector<f64, 2> {
    fn from(p: Point) -> Self {
        Vector([p.x, p.y])
    }
}

#[wasm_bindgen]
pub struct Curve(BezierCurve<f64>);

#[wasm_bindgen]
impl Curve {
    #[wasm_bindgen(constructor)]
    pub fn new(points: Array) -> Curve {
        Curve(BezierCurve(
            (0..points.length()).into_iter()
                .map(|i|
                    points.get(i)
                        .into_serde::<Point>()
                        .unwrap_throw()
                        .into()
                )
                .collect()
        ))
    }

    #[wasm_bindgen]
    pub fn draw(&self, ctx: CanvasRenderingContext2d, color: &JsValue) {
        ctx.begin_path();
        match &self.0[..] {
            &[] | &[_] => {}
            &[a, b] => {
                ctx.move_to(a[0], a[1]);
                ctx.line_to(b[0], b[1]);
            }
            &[a, b, c] => {
                ctx.move_to(a[0], a[1]);
                ctx.quadratic_curve_to(b[0], b[1], c[0], c[1]);
            }
            &[a, b, c, d] => {
                ctx.move_to(a[0], a[1]);
                ctx.bezier_curve_to(b[0], b[1], c[0], c[1], d[0], d[1]);
            }
            &[a, .., b] => {
                ctx.move_to(a[0], a[1]);
                let ticks = 20;
                for i in 1..ticks {
                    let p = self.0.castlejau_eval(i as f64 / ticks as f64);
                    ctx.line_to(p[0], p[1]);
                }
                ctx.line_to(b[0], b[1]);
            }
        }
        ctx.set_stroke_style(color);
        ctx.stroke();
    }

    #[wasm_bindgen(js_name = drawHandles)]
    pub fn draw_handles(&self, ctx: CanvasRenderingContext2d, color: &JsValue) {
        let width = ctx.line_width();
        ctx.set_line_width(width / 2.0);

        for (from, to) in [(self.0[0], self.0[1]), (self.0[self.0.len()-1], self.0[self.0.len()-2])] {
            ctx.begin_path();
            ctx.move_to(from[0], from[1]);
            ctx.line_to(to[0], to[1]);
            ctx.set_stroke_style(color);
            ctx.stroke();
        }
        ctx.set_line_width(width / 4.0);
        ctx.begin_path();
        ctx.move_to(self.0[0][0], self.0[0][1]);
        for p in &self.0[1..] {
            ctx.line_to(p[0], p[1]);
        }
        ctx.stroke();

        ctx.set_line_width(width);
    }

    #[wasm_bindgen(js_name = drawTicks)]
    pub fn draw_ticks(&self, ctx: CanvasRenderingContext2d, color: &JsValue) {
        let width = ctx.line_width();
        for i in 0..=10 {
            let t = i as f64 / 10.0;
            let x = self.0.castlejau_eval(t);
            let dx = self.0.normal(t).normalize() * width;
            let from = x - dx;
            let to = x + dx;
            ctx.begin_path();
            ctx.move_to(from[0], from[1]);
            ctx.line_to(to[0], to[1]);
            ctx.set_stroke_style(color);
            ctx.stroke();
        }
    }

    #[wasm_bindgen(js_name = drawHull)]
    pub fn draw_hull(&self, ctx: CanvasRenderingContext2d, color: &JsValue) {
        let polygon = convex_hull(self.0.iter().map(Clone::clone).collect());

        let width = ctx.line_width();
        ctx.set_line_width(width / 4.0);

        ctx.begin_path();
        ctx.move_to(polygon[0][0], polygon[0][1]);
        for p in &polygon[1..] {
            ctx.line_to(p[0], p[1]);
        }
        ctx.close_path();
        ctx.set_stroke_style(color);
        ctx.stroke();

        ctx.set_line_width(width);
    }

    #[wasm_bindgen(js_name = drawBoundingBox)]
    pub fn draw_bounding_box(&self, ctx: CanvasRenderingContext2d, color: &JsValue) {
        let bb = self.0.bounding_box();

        let width = ctx.line_width();
        ctx.set_line_width(width / 4.0);

        ctx.begin_path();
        ctx.move_to(bb.min[0], bb.min[1]);
        ctx.line_to(bb.max[0], bb.min[1]);
        ctx.line_to(bb.max[0], bb.max[1]);
        ctx.line_to(bb.min[0], bb.max[1]);
        ctx.close_path();
        ctx.set_stroke_style(color);
        ctx.stroke();

        ctx.set_line_width(width);
    }

    #[wasm_bindgen(js_name = drawMinimalBox)]
    pub fn draw_minimal_bounding_box(&self, ctx: CanvasRenderingContext2d, color: &JsValue) {
        if self.0.degree() > 3 {
            return;
        }
        let bb = self.0.minimal_bounding_box();

        let width = ctx.line_width();
        ctx.set_line_width(width / 4.0);

        ctx.begin_path();
        ctx.move_to(bb.min[0], bb.min[1]);
        ctx.line_to(bb.max[0], bb.min[1]);
        ctx.line_to(bb.max[0], bb.max[1]);
        ctx.line_to(bb.min[0], bb.max[1]);
        ctx.close_path();
        ctx.set_stroke_style(color);
        ctx.stroke();

        ctx.set_line_width(width);
    }
}
