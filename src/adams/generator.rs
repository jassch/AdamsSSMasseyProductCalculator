
use super::Bidegree;


use std::fmt;
use std::fmt::{Display, Formatter};

//type AdamsGenerator = (u32, i32, usize);

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct AdamsGenerator {
    /// resolution degree
    s: u32,
    /// internal degree
    t: i32, 
    /// generator index
    idx: usize, 
}

impl AdamsGenerator {
    pub fn s(&self) -> u32 {
        self.s
    }
    pub fn t(&self) -> i32 {
        self.t
    }
    pub fn degree(&self) -> Bidegree {
        (self.s, self.t)
    }
    pub fn n(&self) -> i32 {
        self.t-self.s as i32
    }
    pub fn idx(&self) -> usize {
        self.idx
    }
    pub fn new(s: u32, t: i32, idx: usize) -> AdamsGenerator {
        AdamsGenerator {
            s,
            t,
            idx,
        }
    }
}

impl Display for AdamsGenerator {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        write!(f, "({}, {}, {})", self.t, self.n(), self.idx())
    }
}

impl From<(u32,i32,usize)> for AdamsGenerator {
    fn from(tuple: (u32, i32, usize)) -> Self {
        Self::new(tuple.0, tuple.1, tuple.2)
    }
}
impl From<(Bidegree,usize)> for AdamsGenerator {
    fn from(tuple: (Bidegree, usize)) -> Self {
        let ((s,t), idx) = tuple;
        Self::new(s, t, idx)
    }
}

impl From<AdamsGenerator> for (u32,i32,usize) {
    fn from(gen: AdamsGenerator) -> Self {
        (gen.s(), gen.t(), gen.idx())
    }
}