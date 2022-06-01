use crate::vector::Vector;

pub struct BoundingBox<K: Copy + PartialOrd> {
    pub min: Vector<K, 2>,
    pub max: Vector<K, 2>,
}

impl <T: PartialOrd + Copy> BoundingBox<T> {
    pub fn from_iter<Iter: Iterator<Item=Vector<T, 2>>>(mut points: Iter) -> BoundingBox<T> {
        let mut min = points.next().expect("Should at least contain two point");
        let mut max = min;
        for p in points {
            if min[0] > p[0] {
                min[0] = p[0];
            }
            if min[1] > p[1] {
                min[1] = p[1];
            }
            if max[0] < p[0] {
                max[0] = p [0];
            }
            if max[1] < p[1] {
                max[1] = p[1];
            }
        }
        BoundingBox {min, max}
    }

    pub fn from_slice(points: &[Vector<T, 2>]) -> BoundingBox<T> {
        BoundingBox::from_iter(points.iter().map(|&p| p))
    }

    pub fn contains(&self, point: Vector<T, 2>) -> bool {
        self.min <= point && point <= self.max
    }

    pub fn intersects(&self, other: &Self) -> bool {
        self.intersecting_interval::<0>(other) && self.intersecting_interval::<1>(other)
    }

    fn intersecting_interval<const I: usize>(&self, other: &Self) -> bool {
        (self.min[I] < other.min[I] && other.min[I] < self.max[I])
            || (self.min[I] < other.max[I] && other.max[I] < self.max[I])
            || (other.min[I] < self.min[I] && self.min[I] < other.max[I])
            || (other.min[I] < self.max[I] && self.max[I] < other.max[I])
    }
}

impl <K: Copy + PartialOrd> From<[Vector<K, 2>; 2]> for BoundingBox<K> {
    fn from(array: [Vector<K, 2>; 2]) -> Self {
        BoundingBox { min: array[0], max: array[1] }
    }
}

impl <K: Copy + PartialOrd> From<BoundingBox<K>> for [Vector<K, 2>; 2] {
    fn from(bb: BoundingBox<K>) -> Self {
        [bb.min, bb.max]
    }
}