use super::AdamsElement;
use super::Bidegree;
use crate::affinespace::AffineSpace;
use crate::utils;

//use saveload::{Load, Save};

use std::cmp::{Ordering, PartialOrd};

use std::fmt;
use std::fmt::{Display, Formatter};

use std::io;
use std::io::{Read, Write};

use fp::matrix::Subspace;
use fp::prime::ValidPrime;
use fp::vector::FpVector;

//type AdamsElement = (u32,i32,FpVector);

#[derive(Debug, Clone)]
pub struct MasseyProduct {
    /// resolution degree a.s() + b.s() + c.s() - 1
    s: u32,
    /// internal degree a.t() + b.t() + c.t()
    t: i32,
    /// representative
    rep: FpVector,
    /// indeterminacy data
    l_indet: Subspace,
    r_indet: Subspace,
    indet: Subspace,
}

impl MasseyProduct {
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
        self.t - self.s as i32
    }
    pub fn rep(&self) -> &FpVector {
        &self.rep
    }
    pub fn rep_elt(&self) -> AdamsElement {
        AdamsElement::new(self.s, self.t, self.rep.clone())
    }
    pub fn left_indet(&self) -> &Subspace {
        &self.l_indet
    }
    pub fn right_indet(&self) -> &Subspace {
        &self.r_indet
    }
    pub fn indet(&self) -> &Subspace {
        &self.indet
    }
    pub fn affine(&self) -> AffineSpace {
        AffineSpace::new(self.rep.clone(), self.indet.clone())
    }
    pub fn contains_zero(&self) -> bool {
        return self.indet.contains(self.rep.as_slice());
    }
    pub fn new(
        s: u32,
        t: i32,
        rep: FpVector,
        l_indet: Subspace,
        r_indet: Subspace,
    ) -> MasseyProduct {
        let indet = utils::subspace_sum(&l_indet, &r_indet);
        let rrep = if indet.contains(rep.as_slice()) {
            FpVector::new(rep.prime(), rep.len())
        } else {
            rep
        };
        MasseyProduct {
            s,
            t,
            rep: rrep,
            l_indet,
            r_indet,
            indet,
        }
    }
    pub fn new_ae(rep: &AdamsElement, l_indet: Subspace, r_indet: Subspace) -> MasseyProduct {
        Self::new(rep.s(), rep.t(), rep.vec().clone(), l_indet, r_indet)
    }
}

impl Display for MasseyProduct {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        write!(
            f,
            "({}, {}, {}) + {}",
            self.n(),
            self.s(),
            self.rep(),
            self.indet()
        )
    }
}

impl From<(u32, i32, FpVector, Subspace, Subspace)> for MasseyProduct {
    fn from(tuple: (u32, i32, FpVector, Subspace, Subspace)) -> Self {
        Self::new(tuple.0, tuple.1, tuple.2, tuple.3, tuple.4)
    }
}

impl From<(u32, i32, &FpVector, &Subspace, &Subspace)> for MasseyProduct {
    fn from(tuple: (u32, i32, &FpVector, &Subspace, &Subspace)) -> Self {
        Self::new(
            tuple.0,
            tuple.1,
            tuple.2.clone(),
            tuple.3.clone(),
            tuple.4.clone(),
        )
    }
}

impl From<(&AdamsElement, Subspace, Subspace)> for MasseyProduct {
    fn from(tuple: (&AdamsElement, Subspace, Subspace)) -> Self {
        Self::new_ae(tuple.0, tuple.1, tuple.2)
    }
}
impl From<(&AdamsElement, &Subspace, &Subspace)> for MasseyProduct {
    fn from(tuple: (&AdamsElement, &Subspace, &Subspace)) -> Self {
        Self::new_ae(tuple.0, tuple.1.clone(), tuple.2.clone())
    }
}

impl From<MasseyProduct> for (u32, i32, FpVector, Subspace, Subspace, Subspace) {
    fn from(prod: MasseyProduct) -> Self {
        (
            prod.s,
            prod.t,
            prod.rep,
            prod.l_indet,
            prod.r_indet,
            prod.indet,
        ) // taken by move, so move out
    }
}
/*
impl <'a> From<&'a MasseyProduct> for (u32, i32, &'a FpVector) {
    fn from(elt: &'a AdamsElement) -> Self {
        (elt.s(), elt.t(), elt.vec()) // use method .vec() to avoid moving
    }
}
*/

/*
impl Save for MasseyProduct {
    fn save(&self, buffer: &mut impl Write) -> io::Result<()> {
        self.s.save(buffer)?;
        self.t.save(buffer)?;
        self.rep.save(buffer)?;
        self.l_indet.save(buffer)?;
        self.r_indet.save(buffer)?;
        self.indet.save(buffer)?;
        Ok(())
    }
}

impl Load for MasseyProduct {
    type AuxData = ValidPrime;

    fn load(buffer: &mut impl Read, data: &Self::AuxData) -> io::Result<Self> {
        let s = u32::load(buffer, &())?;
        let t = i32::load(buffer, &())?;
        let rep = FpVector::load(buffer, data)?;
        let l_indet = Subspace::load(buffer, data)?;
        let r_indet = Subspace::load(buffer, data)?;
        let indet = Subspace::load(buffer, data)?;
        Ok(MasseyProduct {
            s,
            t,
            rep,
            l_indet,
            r_indet,
            indet,
        })
    }
}
*/

impl PartialEq<AffineSpace> for MasseyProduct {
    fn eq(&self, other: &AffineSpace) -> bool {
        self.affine() == *other
    }
}
impl PartialEq<MasseyProduct> for AffineSpace {
    fn eq(&self, other: &MasseyProduct) -> bool {
        *self == other.affine()
    }
}
impl PartialOrd<AffineSpace> for MasseyProduct {
    fn partial_cmp(&self, other: &AffineSpace) -> Option<Ordering> {
        self.affine().partial_cmp(other)
    }
}
impl PartialOrd<MasseyProduct> for AffineSpace {
    fn partial_cmp(&self, other: &MasseyProduct) -> Option<Ordering> {
        self.partial_cmp(&other.affine())
    }
}

impl PartialEq for MasseyProduct {
    fn eq(&self, other: &Self) -> bool {
        (self.s == other.s)
            && (self.t == other.t)
            && utils::subspace_equality(&self.l_indet, &other.l_indet)
            && utils::subspace_equality(&self.r_indet, &other.r_indet)
            && self.affine() == other.affine()
    }
}

impl Eq for MasseyProduct {}
