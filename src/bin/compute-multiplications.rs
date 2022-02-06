// cargo run --bin compute-multiplications

use massey::*;

use std::collections::hash_map::HashMap;

use anyhow::{anyhow, Result};
use fp::matrix::Matrix;

use adams::{AdamsGenerator, AdamsMultiplication, Bidegree};

fn callback(
    lhs: AdamsGenerator,
    _max_rhs_deg_computed: Bidegree,
    _matrices: &HashMap<Bidegree, Matrix>,
) -> Result<(), String> {
    println!("Multiplications computed for {}", lhs);
    Ok(())
}

fn main() -> Result<()> {
    let save_file_name = String::from("../massey-prod-calc-data/S_2_resolution.data");
    //let resolution_saves_directory = String::from("../massey-prod-calc-data/S_2_resolution_incremental_data");
    let multiplication_data_directory =
        String::from("../massey-prod-calc-data/S_2_multiplication_data");
    //let massey_product_data_directory = String::from("../massey-prod-calc-data/S_2_massey_prod_data");

    println!("Loading resolution...");
    let mut adams_mult: AdamsMultiplication = AdamsMultiplication::new(
        save_file_name,
        None,
        Some(multiplication_data_directory),
        None,
        None,
    )?;
    let prime = adams_mult.prime();

    //fp::vector::initialize_limb_bit_index_table(prime);

    println!("Computing multiplications...");
    match adams_mult.compute_all_multiplications_callback(true, &mut callback) {
        Ok(_) => {}
        Err(err_info) => {
            eprintln!("{}", err_info);
        }
    }

    Ok(())
}
