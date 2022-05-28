use std::cmp::Ordering;
use num::Num;
use crate::vector::Vector;

enum Turn {
    Left, None, Right,
}

pub type Point<T> = Vector<T, 2>;

fn turn_type<R>(x: Point<R>, y: Point<R>, z: Point<R>) -> Turn
    where R: Copy + Num + PartialOrd
{
    // Compute third component of 3d cross product between xy and xz
    let cross = (y - x).cross(&(z - x));
    match PartialOrd::partial_cmp(&cross, &R::zero())
        .expect("NaN shouldn't happen")
    {
        Ordering::Less => Turn::Right,
        Ordering::Equal => Turn::None,
        Ordering::Greater => Turn::Left,
    }
}

fn sort_angle<R>(origin: Point<R>, points: &mut Vec<Point<R>>)
    where R: Copy + Num + PartialOrd
{
    points.sort_by(|&x, &y| {
        match turn_type(origin, x, y) {
            Turn::Left => Ordering::Less,
            Turn::None => {
                let dx = x - origin;
                let dy = y - origin;
                let dist_x = dx * dx;
                let dist_y = dy * dy;
                PartialOrd::partial_cmp(&dist_x, &dist_y)
                    .expect("NaN shouldn't happen")
            }
            Turn::Right => Ordering::Greater,
        }
    });
}

pub fn convex_hull<R>(mut points: Vec<Point<R>>) -> Vec<Point<R>>
    where R: Copy + Num + PartialOrd
{
    let mut stack = Vec::new();

    // Find point with lowest y-coord (if equal lowest x)
    let mut lowest_x = points[0][0];
    let mut lowest_y = points[0][1];
    let mut lowest_i = 0;
    for (i, p) in points.iter().enumerate().skip(1) {
        if p[1] < lowest_y || (p[1] == lowest_y && p[0] < lowest_x) {
            lowest_x = p[0];
            lowest_y = p[1];
            lowest_i = i;
        }
    }

    // Add lowest point to stack and sort rest by angle
    stack.push(points[lowest_i]);
    points.remove(lowest_i);
    sort_angle(stack[0], &mut points);

    // Populate stack only allowing left turns along the hull's boundary
    for p in points.into_iter() {
        while {
            let l = stack.len();
            l > 1 && matches!(turn_type(stack[l-2], stack[l-1], p), Turn::Right)
        } {
            stack.pop();
        }
        stack.push(p);
    }

    stack
}
