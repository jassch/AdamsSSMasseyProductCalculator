
use super::Bidegree;

use saveload::{Save, Load};

use std::fmt;
use std::fmt::{Display, Formatter};

use std::io::{Read, Write};
use std::io;

use fp::vector::FpVector;


//type AdamsElement = (u32,i32,FpVector);

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct AdamsElement {
    /// resolution degree
    s: u32,
    /// internal degree
    t: i32, 
    /// generator index
    vec: FpVector, 
}

impl AdamsElement {
    pub fn s(&self) -> u32 {
        self.s
    }
    pub fn t(&self) -> i32 {
        self.t
    }
    pub fn degree(&self) -> Bidegree {
        (self.s, self.t).into()
    }
    pub fn n(&self) -> i32 {
        self.t-self.s as i32
    }
    pub fn vec(&self) -> FpVector {
        self.vec.clone()
    }
    pub fn new(s: u32, t: i32, vec: FpVector) -> AdamsElement {
        AdamsElement {
            s,
            t,
            vec,
        }
    }
}

impl Display for AdamsElement {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        write!(f, "({}, {}, {})", self.n(), self.s(), self.vec())
    }
}

impl From<(u32,i32,FpVector)> for AdamsElement {
    fn from(tuple: (u32, i32, FpVector)) -> Self {
        Self::new(tuple.0, tuple.1, tuple.2)
    }
}
impl From<(Bidegree,FpVector)> for AdamsElement {
    fn from(tuple: (Bidegree, FpVector)) -> Self {
        let (deg, idx) = tuple;
        let (s, t) = deg.into();
        Self::new(s, t, idx)
    }
}

impl From<AdamsElement> for (u32,i32,FpVector) {
    fn from(elt: AdamsElement) -> Self {
        (elt.s, elt.t, elt.vec) // taken by move, so move out
    }
}
impl From<&AdamsElement> for (u32,i32,FpVector) {
    fn from(elt: &AdamsElement) -> Self {
        (elt.s(), elt.t(), elt.vec()) // use method .vec() to avoid moving
    }
}

impl Save for AdamsElement {
    fn save(&self, buffer: &mut impl Write) -> io::Result<()> {
        self.s.save(buffer)?;
        self.t.save(buffer)?;
        self.vec.save(buffer)?;
        Ok(())
    }
}

impl Load for AdamsElement {
    type AuxData = <FpVector as Load>::AuxData;

    fn load(buffer: &mut impl Read, data: &Self::AuxData) -> io::Result<Self> {
        let s = u32::load(buffer, &())?;
        let t = i32::load(buffer, &())?;
        let vec = FpVector::load(buffer, data)?;
        Ok(AdamsElement {
            s,
            t,
            vec,
        })
    }
}