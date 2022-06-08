use crate::graham_scan::convex_hull;
use crate::bezier::BezierCurve;
use crate::svg::{SVG, Circle, Line, Path, PathInstructions};

impl SVG {
    pub fn add_bezier(&mut self, curve: &BezierCurve<f64>, mut path: Path) {
        if curve.len() < 2 {
            return;
        }
        let mut push = |instr| path.instructions.push((true, instr));
        push(PathInstructions::MoveTo(curve[0]));
        match curve.len() {
            0 | 1 => unreachable!(),
            2 => push(PathInstructions::LineTo(curve[1])),
            3 => push(PathInstructions::Quadratic(curve[1], curve[2])),
            4 => push(PathInstructions::Cubic(curve[1], curve[2], curve[3])),
            _ => {
                let N = 30;
                for i in 1..N {
                    let t = i as f64 / N as f64;
                    let p = curve.castlejau_eval(t);
                    push(PathInstructions::LineTo(p));
                }
                push(PathInstructions::LineTo(curve[curve.len()-1]))
            },
        }
        self.add_elem(path);
    }

    pub fn debug_bezier(&mut self, curve: &BezierCurve<f64>, color: &'static str) {
        let n = curve.len();

        // Draw Curve itself
        let mut path = Path::default();
        path.stroke_color = color;
        self.add_bezier(curve, path);

        // Draw 11 ticks along the curve
        for i in 0..=10 {
            let t = i as f64 / 10.0;
            let x = curve.castlejau_eval(t);
            let dx = curve.normal(t).normalize();
            let from = x - dx;
            let to = x + dx;
            self.add_elem(Circle { center: from, radius: 0.25, color, });
            self.add_elem(Circle { center: to, radius: 0.25, color, });
            self.add_elem(Line { from, to, width: Some(0.5), color, })
        }

        // Draw control points
        for &p in &curve[1..n-1] {
            self.add_elem(Circle {
                center: p,
                radius: 1.0,
                color,
            });
        }

        // Draw handles for last control points
        self.add_elem(Line {
            from: curve[0],
            to: curve[1],
            width: Some(0.5),
            color,
        });
        self.add_elem(Line {
            from: curve[n-1],
            to: curve[n-2],
            width: Some(0.5),
            color,
        });

        // Draw convex hull
        let polygon = convex_hull(curve.0.iter().map(Clone::clone).collect());

        let mut instructions = Vec::with_capacity(polygon.len() + 1);
        instructions.push((true, PathInstructions::MoveTo(polygon[0].into())));
        for p in polygon.into_iter().skip(1) {
            instructions.push((true, PathInstructions::LineTo(p.into())));
        }
        instructions.push((true, PathInstructions::Close));

        self.add_elem(Path {
            stroke_color: color,
            width: 0.1,
            instructions,
            ..Default::default()
        });

        // Draw bounding box
        let bb = curve.bounding_box();
        let size = bb.max - bb.min;
        self.add_elem(Path {
            stroke_color: color,
            width: 0.1,
            instructions: vec![
                (true, PathInstructions::MoveTo(bb.min)),
                (false, PathInstructions::Horizontal(size[0])),
                (false, PathInstructions::Vertical(size[1])),
                (false, PathInstructions::Horizontal(-size[0])),
                (true, PathInstructions::Close),
            ],
            ..Default::default()
        });
    }
}
