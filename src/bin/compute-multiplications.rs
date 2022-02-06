// cargo run --bin compute-multiplications

use massey::*;

use std::collections::hash_map::HashMap;

use anyhow::Result;
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
    let save_file_name = query::with_default(
        "Save directory",
        "../massey-prod-calc-data/S_2_resolution.data",
        |filename| core::result::Result::<_, std::convert::Infallible>::Ok(String::from(filename)),
    );

    let multiplication_data_directory = query::with_default(
        "Multiplication data directory",
        "../massey-prod-calc-data/S_2_multiplication_data",
        |filename| core::result::Result::<_, std::convert::Infallible>::Ok(String::from(filename)),
    );

    println!("Loading resolution...");
    let mut adams_mult: AdamsMultiplication = AdamsMultiplication::new(
        save_file_name,
        None,
        Some(multiplication_data_directory),
        None,
        None,
    )?;

    println!("Computing multiplications...");
    match adams_mult.compute_all_multiplications_callback(true, &mut callback) {
        Ok(_) => {}
        Err(err_info) => {
            eprintln!("{}", err_info);
        }
    }

    Ok(())
}
