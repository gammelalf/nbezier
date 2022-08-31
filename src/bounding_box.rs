//! A simple geometry primitive

use nalgebra::Vector2;

/// An axis aligned bounding box
///
/// Used as geometry primitive in calculating curves' intersections.
pub struct BoundingBox<T: PartialOrd> {
    /// Corner whose coordinates have the lowest values
    pub min: Vector2<T>,

    /// Corner whose coordinates have the highest values
    pub max: Vector2<T>,
}

impl<T: Clone + PartialOrd> BoundingBox<T> {
    /// Create the smallest bounding box containing all points in an iterator.
    pub fn from_iter<Iter: Iterator<Item = Vector2<T>>>(mut points: Iter) -> BoundingBox<T> {
        let mut min = points.next().expect("Should at least contain two point");
        let mut max = min.clone();
        for p in points {
            if min[0] > p[0] {
                min[0] = p[0].clone();
            }
            if min[1] > p[1] {
                min[1] = p[1].clone();
            }
            if max[0] < p[0] {
                max[0] = p[0].clone();
            }
            if max[1] < p[1] {
                max[1] = p[1].clone();
            }
        }
        BoundingBox { min, max }
    }

    /// Create the smallest bounding box containing all points in an slice.
    pub fn from_slice(points: &[Vector2<T>]) -> BoundingBox<T> {
        BoundingBox::from_iter(points.iter().map(|p| p.clone()))
    }

    /// Check whether a point is contained by the box.
    ///
    /// If the point lies on the boundary, it is said to be contained.
    pub fn contains(&self, point: Vector2<T>) -> bool {
        self.min[0] <= point[0]
            && point[0] <= self.max[0]
            && self.min[1] <= point[1]
            && point[1] <= self.max[1]
    }

    /// Check whether this box intersects with another one.
    pub fn intersects(&self, other: &Self) -> bool {
        self.intersecting_interval::<0>(other) && self.intersecting_interval::<1>(other)
    }

    /// Check if two boxes overlap in their projection on the x or y axis.
    fn intersecting_interval<const I: usize>(&self, other: &Self) -> bool {
        (self.min[I] < other.min[I] && other.min[I] < self.max[I])
            || (self.min[I] < other.max[I] && other.max[I] < self.max[I])
            || (other.min[I] < self.min[I] && self.min[I] < other.max[I])
            || (other.min[I] < self.max[I] && self.max[I] < other.max[I])
    }
}