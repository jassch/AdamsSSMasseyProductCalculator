
use saveload::{Save, Load};

use std::io::{Read, Write};
use std::io;

use std::fmt;
use std::fmt::{Display, Formatter};

use fp::vector::FpVector;
use fp::matrix::Subspace;
use fp::prime::ValidPrime;
use std::cmp::{PartialOrd, Ordering};



#[derive(Debug, Clone)]
pub struct AffineSpace {
    translation: FpVector,
    plane: Subspace
}

impl AffineSpace {
    pub fn new(translation: FpVector, plane: Subspace) -> Self {
        AffineSpace {
            translation,
            plane
        }
    }
    pub fn from_subspace(plane: Subspace) -> Self {
        Self::new(FpVector::new(plane.prime(), plane.ambient_dimension()), plane)
    }
    pub fn from_vector(translation: FpVector) -> Self {
        let sub = Subspace::empty_space(translation.prime(), translation.len());
        Self::new(translation, sub)
    }
    pub fn zero(p: ValidPrime, dim: usize) -> Self {
        Self::new(FpVector::new(p, dim), Subspace::empty_space(p, dim))
    }
    pub fn entire_space(p: ValidPrime, dim: usize) -> Self {
        Self::new(FpVector::new(p, dim), Subspace::entire_space(p, dim))
    }
    pub fn prime(&self) -> ValidPrime {
        self.translation.prime()
    }
    pub fn dimension(&self) -> usize {
        self.plane.dimension()
    }
    pub fn representative(&self) -> &FpVector {
        &self.translation
    }
    pub fn plane(&self) -> &Subspace {
        &self.plane
    }
    pub fn contains(&self, other: &FpVector) -> bool {
        let mut diff = other.clone();
        diff.add(&self.translation, *self.prime()-1);
        self.plane.contains(diff.as_slice())
    }
    pub fn contains_zero(&self) -> bool {
        self.plane.contains(self.translation.as_slice())
    }
}

impl PartialEq for AffineSpace {
    fn eq(&self, other: &Self) -> bool {
        if self.prime() != other.prime() {
            return false;
        }
        if self.dimension() != other.dimension() {
            return false;
        }
        for v in self.plane.basis() {
            if !other.plane.contains(v.as_slice()) {
                return false;
            }
        }
        // now know planes are equal
        self.contains(&other.translation)
    }
}

impl Eq for AffineSpace {}

impl PartialOrd for AffineSpace {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        if self.prime() != other.prime() {
            return None;
        }
        let (smaller,larger,swapped) = if self.dimension() <= other.dimension() {
            (self, other, false)
        } else {
            (other, self, true)
        };
        for v in smaller.plane.basis() {
            if !larger.plane.contains(v.as_slice()) {
                return None; // can't compare
            }
        }
        // now know plane of smaller is subspace of plane of larger
        if larger.contains(&smaller.translation) {
            // know smaller <= larger, == if dimensions are equal
            if self.dimension() == other.dimension() {
                return Some(Ordering::Equal)
            }
            // smaller < larger
            if swapped {
                // smaller is other
                Some(Ordering::Greater)
            } else {
                Some(Ordering::Less)
            }
        } else {
            None
        }

    }
}

impl Display for AffineSpace {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        write!(f, "({} + {})", self.representative(), self.plane())
    }
}
impl Save for AffineSpace {
    fn save(&self, buffer: &mut impl Write) -> io::Result<()> {
        self.translation.save(buffer)?;
        self.plane.save(buffer)?;
        Ok(())
    }
}

impl Load for AffineSpace {
    type AuxData = ValidPrime;

    fn load(buffer: &mut impl Read, data: &Self::AuxData) -> io::Result<Self> {
        let translation = FpVector::load(buffer, data)?;
        let plane = Subspace::load(buffer, data)?;
        Ok(AffineSpace {
            translation,
            plane
        })
    }
}