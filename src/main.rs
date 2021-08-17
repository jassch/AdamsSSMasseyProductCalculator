
use std::cmp::min;
use std::fs::File;
use std::path::Path;
use std::sync::Arc;
use std::clone::Clone;
use std::collections::hash_map::HashMap;
use algebra::module::{Module};
//use error::Error;
use ext::chain_complex::{ChainComplex, FreeChainComplex};
use ext::CCC;
use ext::resolution_homomorphism::ResolutionHomomorphism;
use ext::resolution::Resolution;
use ext::utils::construct;
use fp::prime::ValidPrime;
use fp::matrix::Matrix;
use fp::matrix::Subspace;
use fp::vector::FpVector;
use saveload::Save;

pub mod utils;
use utils::AllVectorsIterator;

/* need to store the products 
 * need to be able to extract Massey productable triples
 * then need to compute the Massey products and store them.
 * should be extensible 
 * note Massey productable triples can involve non generators
 * 
 * Multiplication is a bilinear map 
 * Adams(s1,t1) x Adams(s2,t2) -> Adams(s1+s2,t1+t2)
 * Idea 1: 
 * Store bilinear map as 3d matrix 
 * (normal 2d matrix with entries in Adams(s1+s2,t1+t2))
 * Idea 2:
 * Store as linear map from the tensor product 
 * Adams(s1,t1) \otimes Adams(s2,t2) -> Adams(s1+s2,t1+t2)
 * this is a normal matrix
 * Idea 3: 
 * For each generator x_{s1,t1,i} in (s1, t1) store the matrix 
 * for left multiplication x_{s1,t1,i}
 * (this is what we start with)
 * !!! We'll start with this !!! and adjust as necessary
 * 
 * Goal is to compute pairs 
 * (a,b) \in Adams(s1,t1)\times Adams(s2,t2) such that 
 * mu(a,b) = 0
 * 
 */
type Bidegree = (u32, i32);
type AdamsGenerator = (u32, i32, usize);
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
}

impl AdamsMultiplication {
    pub fn new(res_file_name: String, max_s: u32, max_t: i32) -> error::Result<AdamsMultiplication> {
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
            num_gens.insert((s,t), res.number_of_gens_in_bidegree(s,t));
        }
        
        Ok(AdamsMultiplication {
            resolution: res,
            res_file_name: res_file_name,
            num_gens: num_gens,
            max_s: max_s,
            max_t: max_t,
            multiplication_range_computed:
                HashMap::new(),
            multiplication_matrices: 
                HashMap::new(),
        })
    }

    /// return Arc to the resolution
    pub fn resolution(self: &Self) -> Arc<Resolution<CCC>>{
        self.resolution.clone()
    }

    pub fn prime(self: &Self) -> ValidPrime {
        self.resolution.prime()
    }
    
    /// return copy of the filename
    pub fn res_file_name(self: &Self) -> String {
        self.res_file_name.clone()
    }

    /// return nonmutable reference
    pub fn multiplication_matrices(self: &Self) -> &HashMap<AdamsGenerator, HashMap<Bidegree, Matrix>> {
        &self.multiplication_matrices
    }

    pub fn num_gens(self: &Self, s: u32, t: i32) -> Option<usize> {
        self.num_gens.get(&(s,t)).map(|u| -> usize { u.clone() })
    }
    pub fn num_gens_bidegree(self: &Self, deg: Bidegree) -> Option<usize> {
        let (s,t) = deg;
        self.num_gens(s,t)
    }
    pub fn multiplication_in_bounds(self: &Self, deg1: Bidegree, deg2: Bidegree) -> bool {
        let (s1,t1)=deg1;
        let (s2,t2)=deg2;
        if (t1 < 0) || (t2 < 0) {
            return false;
        } else {
            return (s1+s2 < self.max_s) && (t1 + t2 < self.max_t)
        }
    }

    /// is the bilinear map Adams(deg1) x Adams(deg2) -> Adams(deg1+deg2)
    /// completely computed?
    pub fn multiplication_completely_computed(self: &Self, deg1: Bidegree, deg2: Bidegree) -> bool {
        let (s1,t1) = deg1;
        let (s2,t2) = deg2;
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
            match self.multiplication_range_computed.get(&(s1,t1,index)) {
                Some((s2_max,t2_max)) => {
                    if (s2 > *s2_max) || (t2 > *t2_max) {
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
        
        for i in 0..mult_max_s {
            let module = res.module(i); // ith free module
            for j in 0..mult_max_t {
                let gens = &module.gen_names()[j];
                for (idx,g) in gens.iter().enumerate() {
                    // get the hom for the corresponding
                    // generator
                    // lift map dual to g, F_i -> FF_2
                    let adams_gen: AdamsGenerator = (i, j, idx);
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

                    let mult_range_bidegree: Bidegree = (domain_max_s-i, domain_max_t-j);

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
                    for i2 in 0..mult_with_max_s {
                        //let module2 = res.module(i2); // ith free module
                        for j2 in 0..mult_with_max_t {
                            if res.number_of_gens_in_bidegree(i+i2,j+j2)==0 {
                                continue; // nothing in codomain, multiplication is trivially 0
                            }
                            //let gens2 = &module2.gen_names()[j2];
                            let matrix = hom.get_map(i+i2).hom_k(j2);

                            // convert to fp::matrix::Matrix and store
                            matrix_hashmap.insert((i2,j2), Matrix::from_vec(res.prime(), &matrix));

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

    pub fn left_multiplication_by(self: &Self, l: Bidegree, vec: &FpVector, r: Bidegree) -> Option<Matrix> {
        let p = self.prime();
        let (ls, lt) = l;
        let (rs, rt) = r;
        let (gs, gt) = (ls + rs, lt + rt);
        self.num_gens_bidegree(l).and_then(|dim_l| -> Option<Matrix> {
            assert_eq!(dim_l, vec.len()); // vec has to have same length as number of gens in bidegree l

            self.num_gens_bidegree(r).and_then(|dim_r| -> Option<Matrix> {
                self.num_gens(gs, gt).and_then(|dim_g| -> Option<Matrix> {
                    let mut result = Matrix::new(p, dim_r, dim_g);
                    if (dim_l==0) || (dim_r==0) || (dim_g==0) {
                        return Some(result);
                    } else {
                        for ix in 0..dim_l {
                            let coeff = vec.entry(ix);
                            if coeff == 0 {
                                continue;
                            }
                            let matrix = match self.multiplication_matrices.get(&(ls, lt, ix)) {
                                Some(hm) => {
                                    match hm.get(&r) {
                                        Some(m) => m,
                                        None => {
                                            return None;
                                        }
                                    }
                                },
                                None => {
                                    return None; // couldn't find an important multiplication matrix
                                }
                            };
                            result += /* coeff* */ matrix; // coeff is 1 though, so we're good
                        }
                        Some(result)
                    }
                })
            })
        })
    }

    /// compute all massey products of massey-productable triples (a,b,c)  
    /// all of whose bidegrees are less than max_massey
    pub fn brute_force_compute_all_massey_products(self: &Self, max_massey: Bidegree) {
        let p = self.prime();
        let (max_mass_s, max_mass_t) = max_massey;
        // first identify kernels of left multiplication in this range
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
                            let dim3 = match self.num_gens(s3, t3) {
                                Some(n) => n,
                                None => { continue; } // not computed. this shouldn't happen
                            };
                            if dim2 == 0 {
                                continue; // no nonzero vectors
                            }
                            let lmul_v1 = match self.left_multiplication_by((s1, t1), &v1, (s2, t2)) {
                                Some(m) => m,
                                None => {
                                    continue;
                                }
                            };
                            let (aug_start, mut lmul_v1_aug) = Matrix::augmented_from_vec(p, &lmul_v1.to_vec());
                            lmul_v1_aug.row_reduce();
                            let kernel_lmul_v1 = lmul_v1_aug.compute_kernel(aug_start);
                            for v2 in AllVectorsIterator::new(&kernel_lmul_v1) {
                                if v2.is_zero() {
                                    continue;
                                }
                                println!("({},{},{})*({},{},{}) = ({},{},{})", s1, t1, v1, s2, t2, v2, s3, t3, format!("0_{}", dim3));
                            }
                        }
                    }
                }
            }
        }
    }

}


fn main() -> error::Result {
    let save_file_name = String::from("S_2_resolution.data");
    
    let max_s=30;
    let max_t=60;
    let mult_max_s=15;
    let mult_max_t=30;
    let mult_with_max_s=15;
    let mult_with_max_t=30;

    let mut adams_mult: AdamsMultiplication = AdamsMultiplication::new(save_file_name, max_s, max_t)?;

    fp::vector::initialize_limb_bit_index_table(adams_mult.resolution().prime());

    adams_mult.compute_multiplications(mult_max_s, mult_max_t, mult_with_max_s, mult_with_max_t);
    adams_mult.brute_force_compute_all_massey_products((7,30));

    /*
    println!("Iterate over whole F_2^5");
    for fp_vec in AllVectorsIterator::new_whole_space(adams_mult.prime(), 5) {
        println!("fp_vec: {}", fp_vec);
    }

    let p = adams_mult.prime();
    let input = [ vec![1, 0, 0, 1, 1]
                , vec![0, 1, 0, 1, 0]
                , vec![0, 0, 1, 0, 1]
                ];
    let mut m = Matrix::from_vec(p, &input);
    m.row_reduce();
    let subspace = Subspace {
        matrix: m
    };

    println!("Iterate over subspace");
    for fp_vec in AllVectorsIterator::new(&subspace) {
        println!("fp_vec: {}", fp_vec);
    }
    */


    //adams_mult.possible_nontrivial_massey_products();



    Ok(())
}

#[allow(dead_code)]
fn old_main() -> error::Result {
    let save_path = Path::new("S_2_resolution.data");
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
    //let res_arc = Arc::new(res);
    let max_s=30;
    let max_t=60;
    let mult_max_s=15;
    let mult_max_t=30;
    let mult_with_max_s=15;
    let mult_with_max_t=30;

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
    
    
    //let cx : Arc<CCC> = res.complex();
    for i in 0..max_s {
        let module = res.module(i);
        println!("Module s={}, t=[0,{}]", i, module.max_computed_degree());
        println!("Module gens: {:?}", module.gen_names());
        /*
        for (i, v) in module.gen_names().iter_enum() {
            println!("Module deg {}: {:?}", i, v);
        }
        */
        //for j in 0 .. min(max_t,module.max_computed_degree()) {
        //    println!("Module ({},{}): dimension {}", i, j, module.dimension(j));
        //}
    }

    for i in 0..mult_max_s {
        let module = res.module(i); // ith free module
        for j in 0..mult_max_t {
            let gens = &module.gen_names()[j];
            for (idx,g) in gens.iter().enumerate() {
                // get the hom for the corresponding
                // generator
                // lift map dual to g, F_i -> FF_2
                let hom = 
                    ResolutionHomomorphism::new(
                        format!("mult-by-{}",g),
                        res.clone(),
                        res.clone(),
                        i, // s
                        j  // t
                    );
                // matrix defining the first hom
                let mut matrix = Matrix::new(
                    res.prime(),
                    module.number_of_gens_in_degree(j),
                    1
                );
                matrix[0].set_entry(idx,1);
                // should send the generator to 1
                // and everything else to 0
                hom.extend_step(i, j, Some(&matrix));
                // give it the first map

                let domain_max_s = min(i+mult_with_max_s,max_s);
                let domain_max_t = min(j+mult_with_max_t,max_t);
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
                /*
                #[cfg(not(feature = "concurrent"))]
                hom.extend_all();

                #[cfg(feature = "concurrent")]
                hom.extend_all_concurrent(&bucket);
                */

                // now read off products
                // product of g with g' is
                // given by composing the lift of the ext class
                // dual to g with the ext class dual to 
                // g'
                // and reading off 
                // for now, let's just print based on the lift_hom code

                /*
                println!("hom mult-by-{} of degree ({},{}): ", g, hom.shift_s, hom.shift_t);

                for (s, n, t) in hom.target.iter_stem() {
                    if s + i >= hom.source.next_homological_degree() 
                        || t + j > hom.source.module(s+i).max_computed_degree()
                        || s + i > domain_max_s
                        || t + j > domain_max_t 
                    {
                        
                        continue; // out of range for computed stuff
                    }
                    let matrix = hom.get_map(s+i).hom_k(t);
                    for (i2, r) in matrix.iter().enumerate() {
                        println!("mult-by-{}(x_({}, {}, {})) = {:?}", g, n, s, i2, r);
                    }
                }
                */

                // ok let's do the proper multiplications
                for i2 in 0..mult_with_max_s {
                    let module2 = res.module(i2); // ith free module
                    for j2 in 0..mult_with_max_t {
                        if res.number_of_gens_in_bidegree(i+i2,j+j2)==0 {
                            continue;
                        }
                        let gens2 = &module2.gen_names()[j2];
                        let matrix = hom.get_map(i+i2).hom_k(j2);
                        for (idx2,g2) in gens2.iter().enumerate() {
                            print!("{} in ({},{}) * {} in ({},{}) = ", g, i, j, g2, i2, j2);
                            if matrix[idx2].len() == 0  {
                                println!("0 (trivial)");
                            } else {
                                println!("{:?} in ({},{})", matrix[idx2], i+i2, j+j2);
                            }
                        }
                    }
                }
            }
        }
    }

    Ok(())
}
