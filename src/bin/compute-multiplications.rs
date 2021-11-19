// cargo run --bin compute-multiplications

use massey::*;

use std::cmp::min;
use std::io;
use std::collections::hash_map::HashMap;
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


use adams::{Bidegree, AdamsElement, AdamsGenerator, AdamsMultiplication, MasseyProduct};

fn callback(lhs: AdamsGenerator, max_rhs_deg_computed: Bidegree, matrices: &HashMap<Bidegree, Matrix>) 
    -> Result<(), String>
{
    println!("Multiplications computed for {}", lhs);
    Ok(())
}


fn main() -> error::Result {
    let save_file_name = String::from("../massey-prod-calc-data/S_2_resolution.data");
    //let resolution_saves_directory = String::from("../massey-prod-calc-data/S_2_resolution_incremental_data");
    let multiplication_data_directory = String::from("../massey-prod-calc-data/S_2_multiplication_data");
    //let massey_product_data_directory = String::from("../massey-prod-calc-data/S_2_massey_prod_data");

    
    println!("Loading resolution...");
    let mut adams_mult: AdamsMultiplication = AdamsMultiplication::new(save_file_name, None, Some(multiplication_data_directory), None, None)?;
    let prime = adams_mult.prime();

    fp::vector::initialize_limb_bit_index_table(prime);

    println!("Computing multiplications...");
    match adams_mult.compute_all_multiplications_callback(true, &mut callback) {
        Ok(_) => {},
        Err(err_info) => {
            eprintln!("{}", err_info);
        }
    }

    




    Ok(())
}