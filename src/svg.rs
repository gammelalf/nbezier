//! Small library to render bezier curves as svg.
//!
//! Only used inside examples and not even exposed.

use nalgebra::Vector2;
use std::fmt::{Display, Formatter};
use nbezier::BezierCurve;

type Rect = (f64, f64, f64, f64);

pub struct SVG {
    pub view_box: Rect,
    pub elements: Vec<Box<dyn Display>>,
}

impl SVG {
    pub fn add_elem<E: Display + 'static>(&mut self, elem: E) {
        self.elements.push(Box::new(elem));
    }
}

impl Display for SVG {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        writeln!(
            f,
            "<svg viewBox=\"{} {} {} {}\" xmlns=\"http://www.w3.org/2000/svg\">",
            self.view_box.0, self.view_box.1, self.view_box.2, self.view_box.3
        )?;
        for elem in self.elements.iter() {
            elem.fmt(f)?;
        }
        writeln!(f, "</svg>")?;
        return Ok(());
    }
}

pub struct Line {
    pub from: Vector2<f64>,
    pub to: Vector2<f64>,
    pub width: Option<f64>,
    pub color: &'static str,
}

impl Display for Line {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "<line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" stroke=\"{}\"",
            self.from[0], self.from[1], self.to[0], self.to[1], self.color
        )?;
        if let Some(width) = self.width {
            write!(f, " stroke-width=\"{}\"", width)?;
        }
        writeln!(f, "/>")?;
        return Ok(());
    }
}

pub struct Circle {
    pub center: Vector2<f64>,
    pub radius: f64,
    pub color: &'static str,
}

impl Display for Circle {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "<circle cx=\"{}\" cy=\"{}\" r=\"{}\" fill=\"{}\"",
            self.center[0], self.center[1], self.radius, self.color
        )?;
        writeln!(f, "/>")?;
        return Ok(());
    }
}

pub struct Path {
    pub stroke_color: &'static str,
    pub fill_color: &'static str,
    pub width: f64,
    pub instructions: Vec<(bool, PathInstructions)>,
}
impl Default for Path {
    fn default() -> Self {
        Path {
            stroke_color: "black",
            fill_color: "none",
            width: 1.0,
            instructions: Vec::with_capacity(2),
        }
    }
}

#[allow(unused)]
pub enum PathInstructions {
    MoveTo(Vector2<f64>),
    LineTo(Vector2<f64>),
    Vertical(f64),
    Horizontal(f64),
    Cubic(Vector2<f64>, Vector2<f64>, Vector2<f64>),
    SmoothCubic(Vector2<f64>, Vector2<f64>),
    Quadratic(Vector2<f64>, Vector2<f64>),
    SmoothQuadratic(Vector2<f64>),
    Elliptic {
        radii: Vector2<f64>,
        angle: f64,
        large_arc: bool,
        sweep: bool,
        center: Vector2<f64>,
    },
    Close,
}

impl Display for Path {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "<path stroke=\"{}\" fill=\"{}\" stroke-width=\"{}\" d=\"",
            self.stroke_color, self.fill_color, self.width
        )?;
        for (is_absolute, instruction) in self.instructions.iter() {
            use PathInstructions::*;
            let a = *is_absolute;
            match instruction {
                MoveTo(p) => write!(f, "{} {} {} ", if a { "M" } else { "m" }, p[0], p[1]),
                LineTo(p) => write!(f, "{} {} {} ", if a { "L" } else { "l" }, p[0], p[1]),
                Vertical(y) => write!(f, "{} {} ", if a { "V" } else { "v" }, y),
                Horizontal(x) => write!(f, "{} {} ", if a { "H" } else { "h" }, x),
                Cubic(c1, c2, e) => write!(
                    f,
                    "{} {} {} {} {} {} {} ",
                    if a { "C" } else { "c" },
                    c1[0],
                    c1[1],
                    c2[0],
                    c2[1],
                    e[0],
                    e[1]
                ),
                SmoothCubic(c2, e) => write!(
                    f,
                    "{} {} {} {} {} ",
                    if a { "S" } else { "s" },
                    c2[0],
                    c2[1],
                    e[0],
                    e[1]
                ),
                Quadratic(c, e) => write!(
                    f,
                    "{} {} {} {} {} ",
                    if a { "Q" } else { "q" },
                    c[0],
                    c[1],
                    e[0],
                    e[1]
                ),
                SmoothQuadratic(e) => write!(f, "{} {} {} ", if a { "T" } else { "t" }, e[0], e[1]),
                Elliptic {
                    radii: r,
                    angle: p,
                    large_arc: l,
                    sweep: s,
                    center: c,
                } => write!(
                    f,
                    "{} {} {} {} {} {} {} {} ",
                    if a { "A" } else { "a" },
                    r[0],
                    r[1],
                    p,
                    *l as u8,
                    *s as u8,
                    c[0],
                    c[1]
                ),
                Close => write!(f, "{}", if a { "Z" } else { "z" },),
            }?
        }
        writeln!(f, "\"/>")?;
        Ok(())
    }
}

/* Methods acutally processing bezier curves */
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
                #[allow(non_snake_case)]
                    let N = 30; // tweakable constant
                for i in 1..N {
                    let t = i as f64 / N as f64;
                    let p = curve.castlejau_eval(t);
                    push(PathInstructions::LineTo(p));
                }
                push(PathInstructions::LineTo(curve[curve.len() - 1]))
            }
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
            self.add_elem(Circle {
                center: from,
                radius: 0.25,
                color,
            });
            self.add_elem(Circle {
                center: to,
                radius: 0.25,
                color,
            });
            self.add_elem(Line {
                from,
                to,
                width: Some(0.5),
                color,
            })
        }

        // Draw control points
        for &p in &curve[1..n - 1] {
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
            from: curve[n - 1],
            to: curve[n - 2],
            width: Some(0.5),
            color,
        });

        // Draw convex hull
        let polygon = curve.convex_hull();

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
