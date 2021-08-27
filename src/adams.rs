
mod element;
mod generator;
mod multiplication;

use saveload::{Save, Load};

use std::io::{Read, Write};
use std::io;

use std::fmt;
use std::fmt::{Display, Formatter};

pub use element::AdamsElement;
pub use generator::AdamsGenerator;
pub use multiplication::AdamsMultiplication;

/// type synonym for (s,t) bidegrees

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Bidegree {
    /// resolution degree
    s: u32,
    /// internal degree
    t: i32, 
}

impl Bidegree {
    pub fn s(&self) -> u32 {
        self.s
    }
    pub fn t(&self) -> i32 {
        self.t
    }
    pub fn n(&self) -> i32 {
        self.t-self.s as i32
    }
    pub fn new(s: u32, t: i32) -> Self {
        Self {
            s,
            t,
        }
    }
}

impl Display for Bidegree {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        write!(f, "({}, {})", self.n(), self.s())
    }
}

impl From<(u32,i32)> for Bidegree {
    fn from(tuple: (u32, i32)) -> Self {
        Self::new(tuple.0, tuple.1)
    }
}

impl From<Bidegree> for (u32,i32) {
    fn from(deg: Bidegree) -> Self {
        (deg.s(), deg.t())
    }
}

impl Save for Bidegree {
    fn save(&self, buffer: &mut impl Write) -> io::Result<()> {
        self.s.save(buffer)?;
        self.t.save(buffer)?;
        Ok(())
    }
}

impl Load for Bidegree {
    type AuxData = ();

    fn load(buffer: &mut impl Read, _: &Self::AuxData) -> io::Result<Self> {
        let s = u32::load(buffer, &())?;
        let t = i32::load(buffer, &())?;
        Ok(Bidegree {
            s,
            t,
        })
    }
}