
use std::cmp::min;
use std::io::Write;
use std::fs::{File, DirBuilder};

use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::clone::Clone;
use std::collections::hash_map::HashMap;

use algebra::module::homomorphism::ModuleHomomorphism;
use algebra::module::{Module};

use ext::chain_complex::{ChainComplex, FreeChainComplex, ChainHomotopy};
use ext::CCC;
use ext::resolution_homomorphism::ResolutionHomomorphism;
use ext::resolution::Resolution;
use ext::utils::construct;

use fp::prime::ValidPrime;
use fp::matrix::Matrix;
use fp::matrix::Subspace;
use fp::vector::FpVector;

use std::io;
use saveload::{Save, Load};

use super::{Bidegree, AdamsElement, AdamsGenerator};

use crate::utils::{AllVectorsIterator, LoadHM, SaveHM};

//#[derive(Clone)]
pub struct AdamsMultiplication {
    /// the resolution object
    resolution: Arc<Resolution<CCC>>,
    // an Arc object to copy
    //resolution_arc: Arc<&Resolution<CCC>>,
    /// filename storing resolution data
    res_file_name: String,
    /// convenience hash map to keep track of the 
    /// number of generators in each degree
    num_gens: 
        HashMap<Bidegree,usize>,
    /// max_s w/ dimensions computed
    max_s: u32,
    /// max_t w/ dimensions computed
    max_t: i32,
    /// keeps track of which degrees the multiplication by 
    /// a given basis element (s,t,index) is computed for
    multiplication_range_computed:
        HashMap<AdamsGenerator, Bidegree>,
    /// stores the multiplication matrices for each degree 
    /// where we could compute the multiplication
    multiplication_matrices: 
        HashMap<AdamsGenerator, HashMap<Bidegree, Matrix>>,
    /// directory to store computed multiplications
    /// each file will have the multiplication matrices as computed so far
    /// it's not great. We'll have to see how we can make this more efficient
    /// possibly need to go in to ext and add Save/Load to resolution homomorphisms
    multiplication_data_directory: String,
    known_decompositions:
        HashMap<AdamsElement, Vec<(AdamsElement, AdamsElement)>>
}


impl AdamsMultiplication {

    fn ensure_multiplication_data_directory_exists(&self) -> io::Result<()> {
        DirBuilder::new()
            .recursive(true)
            .create(self.multiplication_data_directory)
    }

    pub fn new(res_file_name: String, multiplication_data_directory: String, max_s: u32, max_t: i32) -> error::Result<AdamsMultiplication> {
        let save_path = Path::new(&res_file_name);
        //let mut res_opt: Result<Resolution<CCC>,Error> = error::from_string("could not construct module");
        let res_opt;
        {
            let prev_save_file = match File::open(save_path) {
                Err(_why) => None,
                Ok(file) => Some(file),
            };
            res_opt = construct("S_2", prev_save_file);
        }
        let res_no_arc : Resolution<CCC> = res_opt?;
        let res = Arc::new(res_no_arc); 
        // Nvm the following, we're going to move into the Arc here now
        // borrow here so we still own the resolution at res_no_arc
        //let res_arc = Arc::new(res);
    
        let save_file: File = File::create(save_path)?;
    
        #[cfg(not(feature = "concurrent"))]
        res.compute_through_bidegree(max_s, max_t);
    
        #[cfg(feature = "concurrent")]
        {
            let bucket = ext::utils::query_bucket();
            res.compute_through_bidegree_concurrent(max_s, max_t, &bucket);
        }
    
        println!("{}", res.graded_dimension_string());
    
        let mut file = std::io::BufWriter::new(save_file);
        res.save(&mut file)?;

        let mut num_gens: HashMap<Bidegree, usize> = HashMap::new();
        for (s,_n,t) in res.iter_stem() {
            num_gens.insert((s,t).into(), res.number_of_gens_in_bidegree(s,t));
        }
        
        let result = AdamsMultiplication {
            resolution: res,
            res_file_name: res_file_name,
            num_gens: num_gens,
            max_s: max_s,
            max_t: max_t,
            multiplication_range_computed:
                HashMap::new(),
            multiplication_matrices: 
                HashMap::new(),
            multiplication_data_directory,
            known_decompositions:
                HashMap::new(),
        };

        result.ensure_multiplication_data_directory_exists()?;

        Ok(result)
    }

    /// return Arc to the resolution
    pub fn resolution(&self) -> Arc<Resolution<CCC>>{
        self.resolution.clone()
    }

    pub fn prime(&self) -> ValidPrime {
        self.resolution.prime()
    }
    
    /// return copy of the filename
    pub fn res_file_name(&self) -> String {
        self.res_file_name.clone()
    }

    /// return nonmutable reference
    pub fn multiplication_matrices(&self) -> &HashMap<AdamsGenerator, HashMap<Bidegree, Matrix>> {
        &self.multiplication_matrices
    }

    pub fn num_gens(&self, s: u32, t: i32) -> Option<usize> {
        self.num_gens.get(&(s,t).into()).map(|u| -> usize { u.clone() })
    }
    pub fn num_gens_bidegree(&self, deg: Bidegree) -> Option<usize> {
        let (s,t) = deg.into();
        self.num_gens(s,t)
    }
    pub fn multiplication_in_bounds(&self, deg1: Bidegree, deg2: Bidegree) -> bool {
        let (s1,t1)=deg1.into();
        let (s2,t2)=deg2.into();
        if (t1 < 0) || (t2 < 0) {
            return false;
        } else {
            return (s1+s2 < self.max_s) && (t1 + t2 < self.max_t)
        }
    }

    pub fn adams_gen_to_resoln_hom(&self, g: AdamsGenerator) -> 
        Result<ResolutionHomomorphism<Resolution<CCC>,Resolution<CCC>>, String> {
        let (s,t,idx) = g.into();
        let hom = ResolutionHomomorphism::new(
            format!("({},{},{})", s, t, idx),
            self.resolution(),
            self.resolution(),
            s,
            t
        );
        let dim = self.num_gens(s, t).ok_or(format!("resolution not computed through ({}, {})", s,t))?;
        let mut matrix = Matrix::new(self.prime(), dim, 1);
        matrix[idx].set_entry(0, 1);
        hom.extend_step(s, t, Some(&matrix));
        Ok(hom)
    }

    pub fn adams_elt_to_resoln_hom(&self, e: &AdamsElement) -> ResolutionHomomorphism<Resolution<CCC>,Resolution<CCC>> {
        let (s,t,v) = e.into();
        let hom = ResolutionHomomorphism::new(
            format!("({},{},{})", s, t, v),
            self.resolution(),
            self.resolution(),
            s,
            t
        );
        let mut matrix = Matrix::new(v.prime(), v.len(), 1);
        for idx in 0..v.len() {
            matrix[idx].set_entry(0,v.entry(idx));
        }
        hom.extend_step(s, t, Some(&matrix));
        hom
    }
    
    fn multiplication_file_name(g: AdamsGenerator) -> String {
        format!("mult_s{}_t{}_{}.data", g.s(), g.t(), g.idx())
    }

    fn multiplication_file_path(&self, g: AdamsGenerator) -> PathBuf {
        let path: PathBuf = 
            [self.multiplication_data_directory, Self::multiplication_file_name(g)]
                .iter()
                .collect();
        path
    }

    pub fn load_multiplications_for(&self, g: AdamsGenerator) 
        -> io::Result<Option<(Bidegree, HashMap<Bidegree, Matrix>)>> {
        let path = self.multiplication_file_path(g);
        if path.exists() {
            let mut file = File::open(path)?;
            let max_s = u32::load(&mut file, &())?;
            let max_t = i32::load(&mut file, &())?;
            let hm: HashMap<Bidegree, Matrix> = LoadHM::load(&mut file, &((), self.prime()))?.into();
            Ok(Some(((max_s, max_t).into(),hm)))
        } else {
            // this is actually fine though
            Ok(None)
        }
        
    }

    pub fn compute_multiplication(&self, g: AdamsGenerator, mult_with_max: Bidegree) 
        -> Result<(Bidegree, HashMap<Bidegree, Matrix>), String> {
        let (s, t, idx) = g.into();
        if s > self.max_s || t > self.max_t {
            return Err(format!("({}, {}, {}) is out of computed range: ({}, {})", s, t, idx, self.max_s, self.max_t));
        }
        let (rmax_s, rmax_t) = mult_with_max.into();
        let rmax_poss_s = self.max_s - s as u32;
        let rmax_poss_t = self.max_t - t as i32;

        let actual_rmax_s = min(rmax_s, rmax_poss_s);
        let actual_rmax_t = min(rmax_t, rmax_poss_t);
        
        let hm = self.load_multiplications_for(g)?;

        let hom = self.adams_gen_to_resoln_hom(g)?;

    }

    /// is the bilinear map Adams(deg1) x Adams(deg2) -> Adams(deg1+deg2)
    /// completely computed?
    pub fn multiplication_completely_computed(self: &Self, deg1: Bidegree, deg2: Bidegree) -> bool {
        let (s1,t1) = deg1.into();
        let (s2,t2) = deg2.into();
        if !self.multiplication_in_bounds(deg1, deg2) {
            return false;
        }
        let num_gens_left = match self.num_gens_bidegree(deg1) {
            Some(n) => n,
            None => { 
                return false; 
            }
        };
        for index in 0..num_gens_left {
            match self.multiplication_range_computed.get(&(s1,t1,index).into()) {
                Some(deg2_max) => {
                    let (s2_max,t2_max) = (*deg2_max).into();
                    if (s2 > s2_max) || (t2 > t2_max) {
                        return false;
                    }
                }
                None => {
                    return false;
                }
            }
        }
        return true;
    }

    pub fn compute_all_multiplications(self: &mut Self) {
        self.compute_multiplications(self.max_s, self.max_t, self.max_s, self.max_t);
    }

    pub fn compute_multiplications(self: &mut Self, mult_max_s:u32, mult_max_t:i32, mult_with_max_s:u32, mult_with_max_t:i32) {
        let mult_max_s = min(mult_max_s, self.max_s);
        let mult_max_t = min(mult_max_t, self.max_t);

        let res=&self.resolution; // convenience alias
        
        for i in 0..mult_max_s+1 {
            let module = res.module(i); // ith free module
            for j in 0..mult_max_t+1 {
                let gens = &module.gen_names()[j];
                for (idx,g) in gens.iter().enumerate() {
                    // get the hom for the corresponding
                    // generator
                    // lift map dual to g, F_i -> FF_2
                    let adams_gen: AdamsGenerator = (i, j, idx).into();
                    let hom = 
                        ResolutionHomomorphism::new(
                            format!("mult-by-{}",g),
                            res.clone(),
                            res.clone(),
                            i, // s
                            j  // t
                        );
                    // matrix defining the first hom
                    let ngens = module.number_of_gens_in_degree(j);
                    let mut matrix = Matrix::new(
                        res.prime(),
                        ngens,
                        1
                    );
                    //assert_eq!(matrix.rows(), ngens);
                    //assert_eq!(matrix.columns(), 1);
                    
                    matrix[idx].set_entry(0,1);
                    // should send the generator to 1
                    // and everything else to 0
                    hom.extend_step(i, j, Some(&matrix));
                    // give it the first map

                    let domain_max_s = min(i+mult_with_max_s,self.max_s);
                    let domain_max_t = min(j+mult_with_max_t,self.max_t);

                    let (mr_s, mr_t) = (domain_max_s-i, domain_max_t-j);
                    let mult_range_bidegree = (mr_s, mr_t).into();

                    // extend hom            
                    #[cfg(not(feature = "concurrent"))]
                    hom.extend(
                        domain_max_s, 
                        domain_max_t
                    );

                    #[cfg(feature = "concurrent")]
                    hom.extend_concurrent(
                        domain_max_s, 
                        domain_max_t, 
                        &bucket);
                    

                    self.multiplication_range_computed.insert(adams_gen, mult_range_bidegree);

                    let mut matrix_hashmap: HashMap<Bidegree, Matrix> = HashMap::new();

                    // ok let's do the proper multiplications
                    for i2 in 0..mr_s+1 {
                        //let module2 = res.module(i2); // ith free module
                        for j2 in 0..mr_t+1 {
                            if !res.has_computed_bidegree(i+i2, j+j2) {
                                continue; 
                            }
                            if res.number_of_gens_in_bidegree(i+i2,j+j2)==0 {
                                continue; // nothing in codomain, multiplication is trivially 0
                            }
                            //let gens2 = &module2.gen_names()[j2];
                            let matrix = hom.get_map(i+i2).hom_k(j2);

                            // convert to fp::matrix::Matrix and store
                            matrix_hashmap.insert((i2,j2).into(), Matrix::from_vec(res.prime(), &matrix));

                            /*
                            for (idx2,g2) in gens2.iter().enumerate() {
                                print!("{} in ({},{}) * {} in ({},{}) = ", g, i, j, g2, i2, j2);
                                if matrix[idx2].len() == 0  {
                                    println!("0 (trivial)");
                                } else {
                                    println!("{:?} in ({},{})", matrix[idx2], i+i2, j+j2);
                                }
                            }
                            */
                        }
                    }
                    self.multiplication_matrices.insert((i,j,idx).into(), matrix_hashmap);
                }
            }
        }
    }

    /// only uses the dimensions, none of the multiplicative structure
    pub fn possible_nontrivial_massey_products(self: &Self) {
        let mut count = 0;
        let mut not_all_dim_1 = 0;
        let mut max_triples = 0;
        for s1 in 1..self.max_s {
            for t1 in s1 as i32..self.max_t {
                let n1 = match self.num_gens(s1, t1) {
                    Some(n) => n,
                    None => { continue; }
                };
                if n1 == 0 {
                    continue; // empty bidegree, massey products will be uninteresting
                }
                for s2 in 1..self.max_s {
                    for t2 in s2 as i32..self.max_t {
                        let n2 = match self.num_gens(s2, t2) {
                            Some(n) => n,
                            None => { continue; }
                        };
                        if n2 == 0 {
                            continue; // empty bidegree, massey products will be uninteresting
                        }

                        let shift_s=s1+s2-1;
                        let shift_t=t1+t2;

                        if shift_s > self.max_s {
                            continue;
                        }

                        for s3 in 1..self.max_s-shift_s {
                            for t3 in s3 as i32..self.max_t-shift_t {
                                let n3 = match self.num_gens(s3, t3) {
                                    Some(n) => n,
                                    None => { continue; }
                                };
                                if n3 == 0 {
                                    continue; // empty bidegree, massey products will be uninteresting
                                }
                                let final_s = shift_s+s3;
                                let final_t = shift_t+t3;
                                let target_n = match self.num_gens(final_s, final_t) {
                                    Some(n) => n,
                                    None => { continue; }
                                };
                                if target_n == 0 {
                                    continue;
                                }
                                if n1*n2*n3 > 1 {
                                    println!("Potential massey products ({},{}) x ({}, {}) x ({}, {}) -> ({}, {})", s1, t1, s2, t2, s3, t3, final_s, final_t);
                                    println!("Dimensions: {} x {} x {} -> {}", n1, n2, n3, target_n);
                                    not_all_dim_1=not_all_dim_1+1;
                                }
                                count = count+1;
                                max_triples = max_triples + n1*n2*n3;
                            }
                        }
                    }
                }
            }
        }
        println!("VS triples: {}", count);
        println!("VS triples not all linear: {}", not_all_dim_1);
        println!("Max elt triples: {}", max_triples);
    }

    pub fn left_multiplication_by(self: &Self, l: Bidegree, vec: &FpVector, r: Bidegree) -> Result<Matrix, String> {
        let p = self.prime();
        let (ls, lt) = l.into();
        let (rs, rt) = r.into();
        let (gs, gt) = (ls + rs, lt + rt);
        self.num_gens_bidegree(l)
            .ok_or(format!("Couldn't get gens in left bidegree ({}, {})", ls , lt))
            .and_then(|dim_l| -> Result<Matrix, String> {
            assert_eq!(dim_l, vec.len()); // vec has to have same length as number of gens in bidegree l

            self.num_gens_bidegree(r)
                .ok_or(format!("Couldn't get gens in right bidegree ({}, {})", rs , rt))
                .and_then(|dim_r| -> Result<Matrix, String> {
                self.num_gens(gs, gt)
                    .ok_or(format!("Couldn't get gens in final bidegree ({}, {})", gs, gt))
                    .and_then(|dim_g| -> Result<Matrix, String> {
                    let mut result = Matrix::new(p, dim_r, dim_g);
                    if (dim_l==0) || (dim_r==0) || (dim_g==0) {
                        return Ok(result);
                    } else {
                        for ix in 0..dim_l {
                            let coeff = vec.entry(ix);
                            if coeff == 0 {
                                continue;
                            }
                            let matrix = match self.multiplication_matrices.get(&(ls, lt, ix).into()) {
                                Some(hm) => {
                                    match hm.get(&r) {
                                        Some(m) => m,
                                        None => {
                                            return Err(format!("Couldn't get multiplication matrix for gen ({}, {}, {}) in bidegree ({}, {})", ls, lt, ix, rs, rt));
                                        }
                                    }
                                },
                                None => {
                                    return Err(format!("Couldn't get multiplication matrices for gen ({}, {}, {})", ls, lt, ix)); // couldn't find an important multiplication matrix
                                }
                            };
                            result += /* coeff* */ matrix; // coeff is 1 though, so we're good
                        }
                        Ok(result)
                    }
                })
            })
        })
    }

    pub fn right_multiplication_by(self: &Self, r: Bidegree, vec: &FpVector, l: Bidegree) -> Result<Matrix, String> {
        // TODO right now I'm assuming that the multiplication in the Adams Spectral Sequence is commutative
        // so we can return
        self.left_multiplication_by(r, vec, l)
        // Check this assumption
    }

    /// Indeterminacy of massey product only depends on the bidegree of the middle term.
    pub fn compute_indeterminacy_of_massey_product(self: &Self, a: &AdamsElement, b: Bidegree, c: &AdamsElement) -> Result<Subspace, String> {
        let (s1, t1, v1) = a.into();
        let (s2, t2) = b.into();
        let (s3, t3, v3) = c.into();
        let (s_left, t_left) = (s1 + s2 -1, t1+t2);
        let (s_right, t_right) = (s2 + s3 -1, t2+t3);
        let (s_tot, t_tot) = (s1 + s2 + s3 - 1, t1 + t2 + t3);
        let p = self.prime();

        let left_dim = match self.num_gens(s_left, t_left) {
            Some(n) => n,
            None => { 
                return Err(format!("Couldn't get dimension of ({}, {})", s_left, t_left));
            }
        };
        let right_dim = match self.num_gens(s_right, t_right) {
            Some(n) => n,
            None => { 
                return Err(format!("Couldn't get dimension of ({}, {})", s_right, t_right));
            }
        };
        let total_dim = match self.num_gens(s_tot, t_tot) {
            Some(n) => n,
            None => { 
                return Err(format!("Couldn't get dimension of ({}, {})", s_tot, t_tot));
            }
        };

        if left_dim == 0 && right_dim == 0 {
            return Ok(Subspace::empty_space(p, total_dim));
        }

        if left_dim == 0 {
            // just compute left multiplication
            
            // from (s_right, t_right) -> (s_tot, t_tot)
            let left_indet = match self.left_multiplication_by((s1, t1).into(), &v1, (s_right, t_right).into()) {
                Ok(lmat) => lmat,
                Err(err) => { 
                    return Err(
                        format!(
                            "Couldn't compute the left multiplication ({}, {}, {})* : ({}, {}) -> ({}, {}) because {}", 
                            s1, t1, v1,
                            s_right, t_right,
                            s_tot, t_tot,
                            err
                        )); 
                }
            };
                
            let (l_aug_start, mut l_indet_aug) = Matrix::augmented_from_vec(p, &left_indet.to_vec());
            l_indet_aug.row_reduce();
                
            return Ok(l_indet_aug.compute_image(left_indet.columns(), l_aug_start));
        }
        
        if right_dim == 0 {
            // just compute left multiplication
            
            // from (s_right, t_right) -> (s_tot, t_tot)
            let right_indet = match self.right_multiplication_by((s3, t3).into(), &v3, (s_left, t_left).into()) {
                Ok(rmat) => rmat,
                Err(err) => { 
                    return Err(
                        format!(
                            "Couldn't compute the right multiplication *({}, {}, {}) : ({}, {}) -> ({}, {}) because {}", 
                            s3, t3, v3,
                            s_left, t_left,
                            s_tot, t_tot,
                            err
                        )); 
                }
            };
                
            let (r_aug_start, mut r_indet_aug) = Matrix::augmented_from_vec(p, &right_indet.to_vec());
            r_indet_aug.row_reduce();
                
            return Ok(r_indet_aug.compute_image(right_indet.columns(), r_aug_start));
        }

        // from (s_right, t_right) -> (s_tot, t_tot)
        let left_indet = match self.left_multiplication_by((s1,t1).into(), &v1, (s_right, t_right).into()) {
            Ok(lmat) => lmat,
            Err(err) => { 
                return Err(
                    format!(
                        "Couldn't compute the left multiplication ({}, {}, {})* : ({}, {}) -> ({}, {}) because {}", 
                        s1, t1, v1,
                        s_right, t_right,
                        s_tot, t_tot,
                        err
                    )); 
            }
        };
        let right_indet = match self.right_multiplication_by((s3, t3).into(), &v3, (s_left, t_left).into()) {
            Ok(rmat) => rmat,
            Err(err) => { 
                return Err(
                    format!(
                        "Couldn't compute the right multiplication *({}, {}, {}) : ({}, {}) -> ({}, {}) because {}", 
                        s3, t3, v3,
                        s_left, t_left,
                        s_tot, t_tot,
                        err
                    )); 
            }
        };
        
        let (l_aug_start, mut l_indet_aug) = Matrix::augmented_from_vec(p, &left_indet.to_vec());
        l_indet_aug.row_reduce();
        let (r_aug_start, mut r_indet_aug) = Matrix::augmented_from_vec(p, &right_indet.to_vec());
        r_indet_aug.row_reduce();
        
        let l_ind_subsp = l_indet_aug.compute_image(left_indet.columns(), l_aug_start);
        let r_ind_subsp = r_indet_aug.compute_image(right_indet.columns(), r_aug_start);
        let mut indet = Subspace::new(p, l_ind_subsp.dimension() + r_ind_subsp.dimension(), l_ind_subsp.ambient_dimension());
        indet.add_vectors(
            l_ind_subsp
                .iter()
                .take(l_ind_subsp.dimension())
                .cloned()
                .chain(
                    r_ind_subsp
                    .iter()
                    .take(r_ind_subsp.dimension())
                    .cloned()));
        Ok(indet)
    }

    /// compute all massey products of massey-productable triples (a,b,c)  
    /// all of whose bidegrees are less than max_massey
    pub fn brute_force_compute_all_massey_products(self: &Self, max_massey: Bidegree) {
        let mut zero_massey_output = match File::create("zero-massey-prods.txt") {
            Err(error) => { eprintln!("Could not open 'zero-massey-prods.txt' for writing: {}", error); return; }
            Ok(file) => file
        };
        let p = self.prime();
        let (max_mass_s, max_mass_t) = max_massey.into();
        // first identify kernels of left multiplication in this range
        let mut kernels: HashMap<(Bidegree,Bidegree), HashMap<FpVector,Subspace>> = HashMap::new();
        for s1 in 1..max_mass_s {
            for t1 in s1 as i32..max_mass_t {
                let dim1 = match self.num_gens(s1, t1) {
                    Some(n) => n,
                    None => { continue; } // not computed. this shouldn't happen
                };
                if dim1 == 0 {
                    continue; // no nonzero vectors
                }
                for v1 in AllVectorsIterator::new_whole_space(p, dim1) {
                    if v1.is_zero() { // no need to consider the 0 vector
                        continue;
                    }
                    // might go out of bounds if max_mass_s, max_mass_t > 0.5 max_s, max_t
                    // TODO
                    for s2 in 1..max_mass_s {
                        for t2 in s2 as i32..max_mass_t {
                            let (s3, t3) = (s1+s2, t1+t2);
                            let dim2 = match self.num_gens(s2, t2) {
                                Some(n) => n,
                                None => { continue; } // not computed. this shouldn't happen
                            };
                            let _dim3 = match self.num_gens(s3, t3) {
                                Some(n) => n,
                                None => { continue; } // not computed. this shouldn't happen
                            };
                            if dim2 == 0 {
                                continue; // no nonzero vectors
                            }
                            let lmul_v1 = match self.left_multiplication_by((s1, t1).into(), &v1, (s2, t2).into()) {
                                Ok(m) => m,
                                Err(_) => {
                                    continue;
                                }
                            };
                            let (aug_start, mut lmul_v1_aug) = Matrix::augmented_from_vec(p, &lmul_v1.to_vec());
                            lmul_v1_aug.row_reduce();
                            let kernel_lmul_v1 = lmul_v1_aug.compute_kernel(aug_start);
                            if kernel_lmul_v1.dimension() == 0 {
                                // kernel trival
                                continue; // skip
                            }
                            let bidegree_pair = ((s1,t1).into(),(s2,t2).into());
                            match kernels.get_mut(&bidegree_pair) {
                                Some(hm) => {
                                    hm.insert(v1.clone(),kernel_lmul_v1.clone());
                                },
                                None => {
                                    let mut new_hm = HashMap::new();
                                    new_hm.insert(v1.clone(), kernel_lmul_v1.clone());
                                    kernels.insert(bidegree_pair, new_hm);
                                }
                            }
                            /*
                            for v2 in AllVectorsIterator::new(&kernel_lmul_v1) {
                                if v2.is_zero() {
                                    continue;
                                }
                                println!("({},{},{})*({},{},{}) = ({},{},{})", s1, t1, v1, s2, t2, v2, s3, t3, format!("0_{}", dim3));
                            }
                            */
                        }
                    }
                }
            }
        }
        // stores interesting (no element zero, and lands in nontrivial degree for now) massey-productable triples
        let mut triples: Vec<(AdamsElement, AdamsElement, AdamsElement)> = Vec::new();
        for s1 in 1..max_mass_s {
            for t1 in s1 as i32..max_mass_t {
                let deg1 = (s1, t1).into();
                for s2 in 1..max_mass_s {
                    for t2 in s2 as i32..max_mass_t {
                        let deg2 = (s2, t2).into();
                        let bideg_pr_l = (deg1, deg2);
                        let hm_kers = match kernels.get(&bideg_pr_l) {
                            Some(hm_kers) => {
                                hm_kers
                            },
                            None => { continue; } // no interesting vectors/kernels in this bidegree pair
                        };
                        for (v1, ker_v1) in hm_kers.iter() {
                            for v2 in AllVectorsIterator::new(&ker_v1) {
                                if v2.is_zero() {
                                    continue; // skip the zero vector
                                }
                                // now iterate over s3, t3
                                for s3 in 1..max_mass_s {
                                    for t3 in s3 as i32..max_mass_t {
                                        let deg3 = (s3, t3).into();
                                        let bideg_pr_r = (deg2, deg3);
                                        let final_bideg = (s1+s2+s3-1, t1+t2+t3).into();
                                        match self.num_gens_bidegree(final_bideg) {
                                            Some(n) => {
                                                // computed bidegree
                                                if n==0 { // but dimension 0, no interesting massey products
                                                    continue;
                                                }
                                            },
                                            None => {
                                                // uncomputed bidegree, skip
                                                continue;
                                            }
                                        }
                                        let hm_kers_2 = match kernels.get(&bideg_pr_r) {
                                            Some(hm_kers_2) => { hm_kers_2 },
                                            None => { continue; } // no interesting vectors/kernels in this 
                                            // bidegree pair
                                        };
                                        let ker_v2 = match hm_kers_2.get(&v2) {
                                            Some(ker_v2) => { ker_v2 },
                                            None => { continue; } // v2 doesn't have an interesting kernel here
                                        };
                                        for v3 in AllVectorsIterator::new(&ker_v2) {
                                            if v3.is_zero() {
                                                continue;
                                            }
                                            triples.push(((s1,t1,v1.clone()).into(), (s2,t2,v2.clone()).into(), (s3, t3, v3.clone()).into()));
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        for (ae1, ae2, ae3) in &triples {
            let (s1,t1,v1) = ae1.into();
            let (s2,t2,v2) = ae2.into();
            let (s3,t3,v3) = ae3.into();
            let (shift_s,shift_t) = (s1+s2-1, t1+t2);
            let shift_n = shift_t-shift_s as i32;
            let (tot_s, tot_t) = (shift_s+s3, shift_t+t3);
            let tot_n = tot_t-tot_s as i32;
            let target_dim = match self.num_gens(tot_s, tot_t) {
                Some(n) => n,
                None => { continue; }
            };

            let res_hom_2 = self.adams_elt_to_resoln_hom(&ae2);
            res_hom_2.extend_through_stem(shift_s, shift_n);
            let res_hom_3 = self.adams_elt_to_resoln_hom(&ae3);
            res_hom_3.extend_through_stem(tot_s, tot_n);

            let homotopy = ChainHomotopy::new(
                &*self.resolution,
                &*self.resolution,
                s2+s3,
                t2+t3,
                |source_s, source_t, idx, row| {
                    let mid_s = source_s - s3;

                    res_hom_3.get_map(source_s)
                        .compose(res_hom_2.get_map(mid_s))
                        .apply_to_basis_element(row.as_slice_mut(), 1, source_t, idx);
                }
            );

            homotopy.extend(tot_s, tot_t);
            let last = homotopy.homotopy(tot_s);
            let mut answer = vec![0; target_dim];
            let mut nonzero = false;
            for i in 0..target_dim {
                let output = last.output(tot_t, i);
                for (k, entry) in v1.iter().enumerate() {
                    if entry != 0 {
                        answer[i] += entry * output.entry(k); // TODO: might need an offset here

                    }
                }
                if answer[i]!=0 { 
                    nonzero=true;
                }
            }
            //for 
            
            if nonzero {
                let massey_rep = FpVector::from_slice(v1.prime(), &answer);
                let indet = match self.compute_indeterminacy_of_massey_product(ae1, (s2, t2).into(), ae3) {
                    Ok(subsp) => subsp,
                    Err(reason) => {
                        println!("< ({}, {}, {}), ({}, {}, {}), ({}, {}, {}) > = ({}, {}, {}) + {:?}", 
                            s1, t1, v1, 
                            s2, t2, v2, 
                            s3, t3, v3,
                            tot_s, tot_t, massey_rep,
                            format!("{} could not compute indeterminacy because {}", "{??}", reason)
                            );
                        // hopefully this doesn't happen
                        continue; // printed out, keep on going
                    }
                };
                print!("< ({}, {}, {}), ({}, {}, {}), ({}, {}, {}) > = ({}, {}, {}) + {:?}", 
                    s1, t1, v1, 
                    s2, t2, v2, 
                    s3, t3, v3,
                    tot_s, tot_t, massey_rep,
                    indet
                );
                if indet.contains(massey_rep.as_slice()) {
                    println!(" = 0 ")
                } else {
                    println!("");
                }
            } else {
                writeln!(zero_massey_output, "< ({}, {}, {}), ({}, {}, {}), ({}, {}, {}) > = 0 + did not compute indeterminacy", 
                    s1, t1, v1, 
                    s2, t2, v2, 
                    s3, t3, v3,
                );
            }
        }
        println!("{} total triples", triples.len());

    }

}
