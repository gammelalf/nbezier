pub mod bezier;

use std::fmt::{Display, Formatter};
use nalgebra::Vector2;

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
        writeln!(f, "<svg viewBox=\"{} {} {} {}\" xmlns=\"http://www.w3.org/2000/svg\">",
               self.view_box.0, self.view_box.1, self.view_box.2, self.view_box.3)?;
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
        write!(f, "<line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" stroke=\"{}\"",
               self.from[0], self.from[1], self.to[0], self.to[1], self.color)?;
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
        write!(f, "<circle cx=\"{}\" cy=\"{}\" r=\"{}\" fill=\"{}\"",
               self.center[0], self.center[1], self.radius, self.color)?;
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
    Elliptic{radii: Vector2<f64>, angle: f64, large_arc: bool, sweep: bool, center: Vector2<f64>},
    Close,
}

impl Display for Path {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "<path stroke=\"{}\" fill=\"{}\" stroke-width=\"{}\" d=\"", self.stroke_color, self.fill_color, self.width)?;
        for (is_absolute, instruction) in self.instructions.iter() {
            use PathInstructions::*;
            let a = *is_absolute;
            match instruction {
                MoveTo(p) =>
                    write!(f, "{} {} {} ", if a {"M"} else {"m"}, p[0], p[1]),
                LineTo(p) =>
                    write!(f, "{} {} {} ", if a {"L"} else {"l"}, p[0], p[1]),
                Vertical(y) =>
                    write!(f, "{} {} ", if a {"V"} else {"v"}, y),
                Horizontal(x) =>
                    write!(f, "{} {} ", if a {"H"} else {"h"}, x),
                Cubic(c1, c2, e) =>
                    write!(f, "{} {} {} {} {} {} {} ", if a {"C"} else {"c"}, c1[0], c1[1], c2[0], c2[1], e[0], e[1]),
                SmoothCubic(c2, e) =>
                    write!(f, "{} {} {} {} {} ", if a {"S"} else {"s"}, c2[0], c2[1], e[0], e[1]),
                Quadratic(c, e) =>
                    write!(f, "{} {} {} {} {} ", if a {"Q"} else {"q"}, c[0], c[1], e[0], e[1]),
                SmoothQuadratic(e) =>
                    write!(f, "{} {} {} ", if a {"T"} else {"t"}, e[0], e[1]),
                Elliptic { radii: r, angle: p, large_arc: l, sweep: s, center: c} =>
                    write!(f, "{} {} {} {} {} {} {} {} ",if a {"A"} else {"a"}, r[0], r[1], p, *l as u8, *s as u8, c[0], c[1]),
                Close =>
                    write!(f, "{}", if a {"Z"} else {"z"}, ),
            }?
        }
        writeln!(f, "\"/>")?;
        Ok(())
    }
}
