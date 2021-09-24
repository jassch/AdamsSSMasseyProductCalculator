
use std::cmp::min;
use std::io;
use std::fs::File;
use std::path::Path;
use std::sync::Arc;
use std::clone::Clone;

use algebra::module::{Module};
//use error::Error;
use ext::chain_complex::{ChainComplex, FreeChainComplex};
use ext::CCC;
use ext::resolution_homomorphism::ResolutionHomomorphism;
use ext::resolution::Resolution;
use ext::utils::construct;
use fp::vector::FpVector;
use fp::matrix::Matrix;
use saveload::Save;

pub mod utils;

pub mod adams;
use adams::{Bidegree, AdamsElement, AdamsGenerator, AdamsMultiplication};

pub mod lattice;

//pub mod computation;
//use computation::ComputationResult;



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



/*
impl Save for AdamsMultiplication {

}
*/




fn main() -> error::Result {
    let save_file_name = String::from("../massey-prod-calc-data/S_2_resolution.data");
    let resolution_saves_directory = String::from("../massey-prod-calc-data/S_2_resolution_incremental_data");
    let multiplication_data_directory = String::from("../massey-prod-calc-data/S_2_multiplication_data");
    let massey_product_data_directory = String::from("../massey-prod-calc-data/S_2_massey_prod_data");

    let max_s=33;
    let max_t=99;
    //let mult_max_s=15;
    //let mult_max_t=30;
    //let mult_with_max_s=15;
    //let mult_with_max_t=30;

    
    println!("Loading and extending resolution...");
    let mut adams_mult: AdamsMultiplication = AdamsMultiplication::new(save_file_name, resolution_saves_directory, multiplication_data_directory, massey_product_data_directory)?;
    let prime = adams_mult.prime();

    adams_mult.extend_resolution_to((max_s,max_t).into())?;

    fp::vector::initialize_limb_bit_index_table(adams_mult.resolution().prime());

    println!("Loading and computing multiplications...");
    match adams_mult.compute_all_multiplications() {
        Ok(_) => {},
        Err(err_info) => {
            eprintln!("{}", err_info);
        }
    }

    let h0 = (1,1,FpVector::from_slice(prime, &vec![1])).into();
    let h1 = (1,2,FpVector::from_slice(prime, &vec![1])).into();
    let max_massey_deg = (32,96).into(); //(25,60).into();//(32, 96).into();
    // compute Massey products 
    // <-,h0,h1>
    println!("Computing kernels for multiplication by h0 = {}...", h0);
    // first compute kernels for h0
    let kers_h0 = match adams_mult.compute_kernels_right_multiplication(&h0, max_massey_deg) {
        Ok(kers) => kers,
        Err(err_info) => {
            eprintln!("{}", err_info);
            // fail
            return error::from_string(err_info);
            //std::process::exit(-1);
        }
    };
    println!("Computing massey products <-,{},{}>...", h0, h1);
    let (deg_computed, massey_h0_h1) = adams_mult.compute_massey_prods_for_pair(&kers_h0, max_massey_deg, &h0, &h1);
    println!("Massey products <-,{},{}> computed through degree {} out of {}", h0, h1, deg_computed, max_massey_deg);
    let shift_deg = (1,3).into();
    for (a, rep, indet) in massey_h0_h1 {
        let rep_ae: AdamsElement = (a.degree() + shift_deg, rep).into();
        println!("<{}, h0, h1> = {} + {}", a, rep_ae, indet.matrix);
    }
    

    //adams_mult.compute_multiplications(mult_max_s, mult_max_t, mult_with_max_s, mult_with_max_t);
    //adams_mult.brute_force_compute_all_massey_products((7,30).into());

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
