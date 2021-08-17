
use fp::prime::ValidPrime;
//use fp::matrix::Matrix;
use fp::matrix::Subspace;
use fp::vector::FpVector;


/// converts a vector in subspace coordinates to global coordinates
pub fn subspace_to_global(subspace: &Subspace, vec: &FpVector) -> FpVector {
    let mut result = FpVector::new(subspace.prime(),subspace.ambient_dimension());
    subspace.apply(result.as_slice_mut(), 1, vec.as_slice());
    result
}


#[derive(Clone, Debug)]
pub struct AllVectorsIterator {
    subspace: Subspace,
    initial: bool,
    current: FpVector,
    start: FpVector,
}

impl AllVectorsIterator {
    pub fn new_from(subspace: &Subspace, start: &FpVector) -> AllVectorsIterator {
        AllVectorsIterator {
            subspace: subspace.clone(),
            initial: true,
            current: start.clone(),
            start: start.clone(),
        }
    }
    pub fn new_whole_space_from(start: &FpVector) -> AllVectorsIterator {
        Self::new_from(&Subspace::entire_space(start.prime(), start.len()), start)
    }
    pub fn new(subspace: &Subspace) -> AllVectorsIterator {
        Self::new_from(subspace, &FpVector::new(subspace.prime(), subspace.dimension()))
    }
    pub fn new_whole_space(p: ValidPrime, dimension: usize) -> AllVectorsIterator {
        Self::new(&Subspace::entire_space(p, dimension))
    }
}

impl Iterator for AllVectorsIterator {
    type Item = FpVector;
    fn next(&mut self) -> Option<Self::Item> {
        if self.initial {
            self.initial=false;
        } else {
            for ix in 0..self.subspace.dimension() {
                if self.current.entry(ix) != *self.subspace.prime()-1  {
                    self.current.set_entry(ix,self.current.entry(ix)+1);
                    break;
                } else {
                    self.current.set_entry(ix,0);
                }
            }
            // advance current
            if self.current == self.start {
                // resets
                self.initial=true;
                return None;
            }
        }
        return Some(subspace_to_global(&self.subspace,&self.current));
    }
}




