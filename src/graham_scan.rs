//! Implementation of [Graham scan](https://en.wikipedia.org/wiki/Graham_scan) for constructing
//! convex hulls.
//!
//! The convex hull of any geometric shape is the smallest possible convex set containing the
//! entire shape. For a set of points this is a polygon.

use nalgebra::{RealField, Vector2, Vector3};
use std::cmp::Ordering;

/// Different types of turns
pub enum Turn {
    Left,
    None,
    Right,
}

/// Identifies the type of turn 3 points form.
pub fn turn_type<F: RealField>(x: &Vector2<F>, y: &Vector2<F>, z: &Vector2<F>) -> Turn {
    // Compute third component of 3d cross product between xy and xz
    let x = Vector3::new(x[0].clone(), x[1].clone(), F::zero());
    let y = Vector3::new(y[0].clone(), y[1].clone(), F::zero());
    let z = Vector3::new(z[0].clone(), z[1].clone(), F::zero());
    let cross = (&y - &x).cross(&(z - x));

    match PartialOrd::partial_cmp(&cross[2], &F::zero()).expect("NaN shouldn't happen") {
        Ordering::Less => Turn::Right,
        Ordering::Equal => Turn::None,
        Ordering::Greater => Turn::Left,
    }
}

/// Sorts a vector of points in increasing order of the angle they and the point `origin` make with the x-axis.
fn sort_angle<F: RealField>(origin: &Vector2<F>, points: &mut Vec<Vector2<F>>) {
    points.sort_by(|x, y| match turn_type(origin, x, y) {
        Turn::Left => Ordering::Less,
        Turn::None => {
            let dx = x - origin;
            let dy = y - origin;
            let dist_x = dx.dot(&dx);
            let dist_y = dy.dot(&dy);
            PartialOrd::partial_cmp(&dist_x, &dist_y).expect("NaN shouldn't happen")
        }
        Turn::Right => Ordering::Greater,
    });
}

/// Computes the convex hull, a polygon, for a set of points.
///
/// The polygon is given as a set of its vertecies in counterclockwise order starting at the lowest
/// one.
pub fn convex_hull<F: RealField>(mut points: Vec<Vector2<F>>) -> Vec<Vector2<F>> {
    let mut stack = Vec::new();

    // Find point with lowest y-coord (if equal lowest x)
    let mut lowest_x = points[0][0].clone();
    let mut lowest_y = points[0][1].clone();
    let mut lowest_i = 0;
    for (i, p) in points.iter().enumerate().skip(1) {
        if p[1] < lowest_y || (p[1] == lowest_y && p[0] < lowest_x) {
            lowest_x = p[0].clone();
            lowest_y = p[1].clone();
            lowest_i = i;
        }
    }

    // Add lowest point to stack and sort rest by angle
    stack.push(points[lowest_i].clone());
    points.remove(lowest_i);
    sort_angle(&stack[0], &mut points);

    // Populate stack only allowing left turns along the hull's boundary
    for p in points.into_iter() {
        while {
            let l = stack.len();
            l > 1 && matches!(turn_type(&stack[l - 2], &stack[l - 1], &p), Turn::Right)
        } {
            stack.pop();
        }
        stack.push(p);
    }

    stack
}
