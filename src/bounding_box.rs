use nalgebra::Vector2;

pub struct BoundingBox<K: PartialOrd> {
    pub min: Vector2<K>,
    pub max: Vector2<K>,
}

impl <T: Clone + PartialOrd> BoundingBox<T> {
    pub fn from_iter<Iter: Iterator<Item=Vector2<T>>>(mut points: Iter) -> BoundingBox<T> {
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
        BoundingBox {min, max}
    }

    pub fn from_slice(points: &[Vector2<T>]) -> BoundingBox<T> {
        BoundingBox::from_iter(points.iter().map(|p| p.clone()))
    }

    pub fn contains(&self, point: Vector2<T>) -> bool {
        self.min[0] <= point[0] && point[0] <= self.max[0]
            && self.min[1] <= point[1] && point[1] <= self.max[1]
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

impl <K: PartialOrd> From<[Vector2<K>; 2]> for BoundingBox<K> {
    fn from(array: [Vector2<K>; 2]) -> Self {
        let [min, max] = array;
        BoundingBox { min, max }
    }
}

impl <K: PartialOrd> From<BoundingBox<K>> for [Vector2<K>; 2] {
    fn from(bb: BoundingBox<K>) -> Self {
        [bb.min, bb.max]
    }
}