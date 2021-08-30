
/// PartialOrder with greatest lower bounds, or meets
pub trait MeetSemilattice: PartialOrd<Self> {
    fn meet(self, rhs: Self) -> Self;
}

pub fn meet<T>(v1: T, v2: T) -> T
    where T: MeetSemilattice {
    T::meet(v1, v2)
}

/// Partial Order with least upper bounds, or joins
pub trait JoinSemilattice: PartialOrd<Self> {
    fn join(self, rhs: Self) -> Self;
}

pub fn join<T>(v1: T, v2: T) -> T
    where T: JoinSemilattice {
    T::join(v1, v2)
}

/// A lattice is simultaneously a MeetSemilattice and a JoinSemilattice
pub trait Lattice: MeetSemilattice + JoinSemilattice {
}

impl <T: Ord> MeetSemilattice for T {
    fn meet(self, rhs: Self) -> Self {
        Self::min(self, rhs)
    }
}
impl <T: Ord> JoinSemilattice for T {
    fn join(self, rhs: Self) -> Self {
        Self::max(self, rhs)
    }
}
impl <T: Ord> Lattice for T {
}