
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
