
use ext::chain_complex::{ChainComplex, FreeChainComplex};
use ext::utils::construct;
use saveload::Save;
use std::fs::File;
use std::path::Path;
use error::Error;
use ext::resolution::Resolution;
use ext::CCC;
use std::sync::Arc;
use std::cmp::min;
use algebra::module::Module;

fn main() -> error::Result {
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
    let res : Resolution<CCC> = res_opt?;
    let max_s=30;
    let max_t=60;

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
        println!("Module: {}", module.to_minimal_json());
        for j in 0 .. min(max_t,module.max_computed_degree()) {
            println!("Module ({},{}): dimension {}", i, j, module.dimension(j));
        }
    }

    Ok(())
}
