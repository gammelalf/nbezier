//! Implementation of [Graham scan](https://en.wikipedia.org/wiki/Graham_scan) for constructing
//! convex hulls.
//!
//! The convex hull of any geometric shape is the smallest possible convex set containing the
//! entire shape. For a set of points this is a polygon.

use nalgebra::{RealField, Vector2, Vector3};
use std::cmp::Ordering;

/// Different types of turns
///
/// This enum is returned by [`turn_type`] as a more readable duplicate of [`Ordering`].
pub enum Turn {
    /// A left turn
    ///
    /// i.e. cross product > 0
    Left,

    /// No turn
    ///
    /// i.e. cross product = 0
    None,

    /// A right turn
    ///
    /// i.e. cross product < 0
    Right,
}

/// Identifies the turn 3 points form by computing the cross product of their differences.
pub fn turn_type<T: RealField>(x: &Vector2<T>, y: &Vector2<T>, z: &Vector2<T>) -> Turn {
    // Compute third component of 3d cross product between xy and xz
    let x = Vector3::new(x[0].clone(), x[1].clone(), T::zero());
    let y = Vector3::new(y[0].clone(), y[1].clone(), T::zero());
    let z = Vector3::new(z[0].clone(), z[1].clone(), T::zero());
    let cross = (&y - &x).cross(&(z - x));

    match PartialOrd::partial_cmp(&cross[2], &T::zero()).expect("NaN shouldn't happen") {
        Ordering::Less => Turn::Right,
        Ordering::Equal => Turn::None,
        Ordering::Greater => Turn::Left,
    }
}

/// Sorts a vector of points in increasing order of the angle they and the point `origin` make with the x-axis.
fn sort_angle<T: RealField>(origin: &Vector2<T>, points: &mut Vec<Vector2<T>>) {
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
pub fn convex_hull<T: RealField>(mut points: Vec<Vector2<T>>) -> Vec<Vector2<T>> {
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
