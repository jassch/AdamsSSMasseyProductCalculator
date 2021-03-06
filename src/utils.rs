use fp::prime::ValidPrime;
//use fp::matrix::Matrix;
use fp::matrix::Subspace;
use fp::vector::FpVector;

//use saveload::{Load, Save};

use std::hash::Hash;
use std::io;
use std::io::{Read, Write};
use std::sync::Arc;

use ext::chain_complex::ChainComplex;
use ext::resolution::Resolution;
use ext::CCC;

use std::collections::HashMap;

/// converts a vector in subspace coordinates to global coordinates
pub fn subspace_to_global(subspace: &Subspace, vec: &FpVector) -> FpVector {
    let mut result = FpVector::new(subspace.prime(), subspace.ambient_dimension());
    subspace.apply(result.as_slice_mut(), 1, vec.as_slice());
    result
}

pub fn subspace_sum(lsub: &Subspace, rsub: &Subspace) -> Subspace {
    let mut indet = Subspace::new(
        lsub.prime(),
        lsub.dimension() + rsub.dimension(),
        lsub.ambient_dimension(),
    );
    indet.add_vectors(
        lsub.iter()
            .take(lsub.dimension())
            .cloned()
            .chain(rsub.iter().take(rsub.dimension()).cloned()),
    );
    indet
}

pub fn subspace_equality(lsub: &Subspace, rsub: &Subspace) -> bool {
    if lsub.prime() != rsub.prime() {
        return false;
    }
    if lsub.dimension() != rsub.dimension() {
        return false;
    }
    for v in lsub.basis() {
        if !rsub.contains(v.as_slice()) {
            return false;
        }
    }
    true
}

pub fn get_max_defined_degree(res: Arc<Resolution<CCC>>) -> (u32, i32) {
    let mut s = 0;
    let mut t = 0;
    while res.has_computed_bidegree(s + 1, t) {
        s += 1;
    }
    while res.has_computed_bidegree(s, t + 1) {
        t += 1;
    }
    (s, t)
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
        Self::new_from(
            subspace,
            &FpVector::new(subspace.prime(), subspace.dimension()),
        )
    }
    pub fn new_whole_space(p: ValidPrime, dimension: usize) -> AllVectorsIterator {
        Self::new(&Subspace::entire_space(p, dimension))
    }
}

impl Iterator for AllVectorsIterator {
    type Item = FpVector;
    fn next(&mut self) -> Option<Self::Item> {
        if self.initial {
            self.initial = false;
        } else {
            for ix in 0..self.subspace.dimension() {
                if self.current.entry(ix) != *self.subspace.prime() - 1 {
                    self.current.set_entry(ix, self.current.entry(ix) + 1);
                    break;
                } else {
                    self.current.set_entry(ix, 0);
                }
            }
            // advance current
            if self.current == self.start {
                // resets
                self.initial = true;
                return None;
            }
        }
        Some(subspace_to_global(&self.subspace, &self.current))
    }
}

/// save takes a reference to a data structure to store
pub struct SaveHM<'a, K, V>(pub &'a HashMap<K, V>);

impl<'a, K, V> From<&'a HashMap<K, V>> for SaveHM<'a, K, V> {
    fn from(hm: &'a HashMap<K, V>) -> Self {
        SaveHM(hm)
    }
}

impl<'a, K, V> From<SaveHM<'a, K, V>> for &'a HashMap<K, V> {
    fn from(slhm: SaveHM<'a, K, V>) -> Self {
        let SaveHM(hm) = slhm;
        hm
    }
}

/*
impl<'a, K: Save, V: Save> Save for SaveHM<'a, K, V> {
    fn save(&self, buffer: &mut impl Write) -> io::Result<()> {
        let SaveHM(hm) = self;
        hm.len().save(buffer)?;
        for (k, v) in hm.iter() {
            k.save(buffer)?;
            v.save(buffer)?;
        }
        Ok(())
    }
}
*/

/// Load returns a new owned data structure
pub struct LoadHM<K, V>(pub HashMap<K, V>);

impl<K, V> From<HashMap<K, V>> for LoadHM<K, V> {
    fn from(hm: HashMap<K, V>) -> Self {
        LoadHM(hm)
    }
}

impl<K, V> From<LoadHM<K, V>> for HashMap<K, V> {
    fn from(slhm: LoadHM<K, V>) -> Self {
        let LoadHM(hm) = slhm;
        hm
    }
}

/*
impl<K: Load + Eq + Hash, V: Load> Load for LoadHM<K, V> {
    type AuxData = (K::AuxData, V::AuxData);
    fn load(buffer: &mut impl Read, data: &Self::AuxData) -> io::Result<Self> {
        let len = usize::load(buffer, &())?;

        let mut result: HashMap<K, V> = HashMap::new();
        for _idx in 0..len {
            let k = K::load(buffer, &(*data).0)?;
            let v = V::load(buffer, &(*data).1)?;
            result.insert(k, v);
        }
        Ok(LoadHM(result))
    }
}
*/
