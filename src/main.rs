
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
use algebra::module::{BoundedModule, Module};
use ext::resolution_homomorphism::ResolutionHomomorphism;
use fp::matrix::Matrix;


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
    let res_no_arc : Resolution<CCC> = res_opt?;
    let res = Arc::new(res_no_arc);
    //let res_arc = Arc::new(res);
    let max_s=30;
    let max_t=60;
    let mult_max_s=10;
    let mult_max_t=10;

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

                // extend hom            
                #[cfg(not(feature = "concurrent"))]
                hom.extend_all();

                #[cfg(feature = "concurrent")]
                hom.extend_all_concurrent(&bucket);

                // now read off products
                // products are given by figuring
                // out where the multiplication
                // map 
            }
        }
    }

    Ok(())
}
