//! Circuit cost model.
//! usage: ${binary_file} estimate ${k}
use std::{
    time::Instant, io, fs, collections::BTreeMap,
};

use crate::{
    arithmetic::{Field, CurveAffine, Engine},
    circuit::{Cell, Layouter, SimpleFloorPlanner},
    plonk::*,
    poly::{commitment::Params, commitment::ParamsVerifier, EvaluationDomain, Rotation},
    transcript::{Blake2bRead, Blake2bWrite, Challenge255},
};
use group::{prime::PrimeCurveAffine, GroupEncoding};
use pairing::bn256::{Bn256, Fr as Fp, G1Affine};
use rand_core::OsRng;

fn measure_elapsed_time<T, F: FnOnce() -> T>(f: F) -> (f64, T) {
    let start = Instant::now();
    let res = f();
    (start.elapsed().as_secs_f64(), res)
}

/// EstimateResult is to store the output of estimate()
#[derive(Debug)]
pub struct EstimateResult {
    prover_time: f64,
    mem_usage: f64,
    aggregate_circuit_size: usize,
}

impl EstimateResult {
    /// print estimation result.
    pub fn print(&self) {
        println!("prover time = {} (s)", self.prover_time);
        println!("memory usage = {} (KB)", self.mem_usage);
        // println!("aggregate circuit size = {}", aggregate_circuit_size);
    }
}

/// estimate is to estimate the prover time, peek memory usage and aggregate circuit size.
pub fn estimate<
    E: Engine,
    ConcreteCircuit: Circuit<E::Scalar>,
>(
    circuit: ConcreteCircuit, 
    res_1: SimLinearResult,
    res_2: SimLinearResult,
    k: usize,
) -> EstimateResult {
    // NOTE(sphere): init params
    // Initialize the polynomial commitment parameters
    let mut cs = ConstraintSystem::default();
    let config = ConcreteCircuit::configure(&mut cs);

    let generate_fake_params = |k| {
        let s = E::Scalar::random(OsRng);
        let rand_c1 = <E::G1Affine as PrimeCurveAffine>::generator() * s;
        let rand_c2 = <E::G2Affine as PrimeCurveAffine>::generator() * s;
        let rand_c1: E::G1Affine = rand_c1.into();
        let n = 1 << k;
        Params {
            k: k as u32,
            n: n as u64,
            g: (0..n).map(|_| rand_c1.clone()).collect(),
            g_lagrange: (0..n).map(|_| rand_c1.clone()).collect(),
            additional_data: Vec::from(rand_c2.to_bytes().as_ref()),
        }
    };

    let estimate_pt_non_linear = |k| {

        let params = generate_fake_params(k);
        
        // Initialize the domain
        let domain = EvaluationDomain::fake_new(cs.degree() as u32, params.k, E::Scalar::random(OsRng));
        
        // NOTE(sphere): count function call
        let FuncCount { num_fft, num_extended_fft, num_msm, num_btree } = dummy_proof::<E::G1Affine, ConcreteCircuit>(
            &params,
            &cs,
            &domain,
        );
        
        let n = 1 << k as usize;
        let rand_vec: Vec::<E::Scalar> = (0..n).map(|_| E::Scalar::random(&mut OsRng)).collect();
        let rand_vec2 = rand_vec.clone();
        let rand_values = domain.lagrange_from_vec(rand_vec);

        // NOTE(sphere): estimate opr time
        //      msm
        let (time_msm, _) = measure_elapsed_time(|| params.commit_lagrange(&rand_values));
        //      fft
        let (time_fft, rand_poly) = measure_elapsed_time(|| domain.lagrange_to_coeff(rand_values));
        //      extended fft
        let (time_extended_fft, _) = measure_elapsed_time(|| domain.coeff_to_extended(rand_poly));
        //      BTree time cost in lookup argument
        let (time_btree, _) = measure_elapsed_time(|| {
            let mut leftover_table_map: BTreeMap<E::Scalar, u32> = rand_vec2
            .iter().take(n)
            .fold(BTreeMap::new(), |mut acc, coeff| {
                *acc.entry(*coeff).or_insert(0) += 1;
                acc
            });
            for item in rand_vec2 {
                if let Some(count) = leftover_table_map.get_mut(&item) {
                    *count -= 1;
                }
            }
        });
        println!("num_fft = {}, time_fft = {}", num_fft, time_fft);
        println!("num_extended_fft = {}, time_extended_fft = {}", num_extended_fft, time_extended_fft);
        println!("num_msm = {}, time_msm = {}", num_msm, time_msm);
        println!("num_btree = {}, time_btree = {}", num_fft, time_fft);
        
        let pt_non_linear = (num_fft as f64) * time_fft +
                             (num_extended_fft as f64) * time_extended_fft +
                             (num_msm as f64) * time_msm +
                             (num_btree as f64) * time_btree;
        println!("pt_non_linear = {}", pt_non_linear);
        println!("");
        pt_non_linear
    };

    let prover_time = estimate_pt_non_linear(k);

    let calc_linear_term = |x_1: f64, y_1: f64, x_2: f64, y_2: f64, x_3 :f64| {
        y_1 + (y_2 - y_1) / (x_2 - x_1) * (x_3 - x_1)
    };

    let mem_usage = calc_linear_term(
        (1 << res_1.k) as f64, res_1.mem_usage,
        (1 << res_2.k) as f64, res_2.mem_usage,
        (1 << k) as f64,
    );

    // NOTE(sphere): calculate aggregate_circuit_size
    
    EstimateResult {
        prover_time,
        mem_usage,
        aggregate_circuit_size: 0,
    }
}

/// SimLinearResult is to store the result of simulate.
#[derive(Debug)]
pub struct SimLinearResult {
    k: usize,
    // prover_time: f64,
    mem_usage: f64,
}

impl SimLinearResult {
    /// read is to read SimLinearResult from a file.
    pub fn read(filepath: String) -> SimLinearResult {
        let data = fs::read_to_string(filepath).expect("read failed");
        let mut data = data.split_whitespace();
        let k = data.next().unwrap().parse().expect("k parse failed");
        // let prover_time = data.next().unwrap().parse().expect("prover_time (s) parse failed");
        let mem_usage = data.next().unwrap().parse().expect("mem_usage (KB) parse in failed");
        SimLinearResult {
            k,
            // prover_time,
            mem_usage,
        }
    }

    /// create new `SimLinearResult`.
    pub fn new(k: usize, mem_usage: f64) -> Self {
        SimLinearResult {
            k,
            mem_usage,
        }
    }
}

/// simulate_circuit is to run a circuit proving process.
pub fn simulate_circuit<
    E: Engine,
    ConcreteCircuit: Circuit<E::Scalar>,
>(circuit: ConcreteCircuit, k: usize) {
    // let public_inputs_size = 0;

    // Initialize the polynomial commitment parameters
    let params: Params<E::G1Affine> = Params::<E::G1Affine>::unsafe_setup::<E>(k as u32);
    // let params_verifier: ParamsVerifier<E> = params.verifier(public_inputs_size).unwrap();

    // Initialize the proving key
    let vk = keygen_vk(&params, &circuit).expect("keygen_vk should not fail");
    let pk = keygen_pk(&params, vk, &circuit).expect("keygen_pk should not fail");

    // Create a proof
    let mut transcript = Blake2bWrite::<_, _, Challenge255<_>>::init(vec![]);

    let (prover_time, _) = measure_elapsed_time(|| 
        create_proof(&params, &pk, &[circuit], &[&[]], OsRng, &mut transcript)
        .expect("proof generation should not fail")
    );

    // NOTE(liutainyi): output prover_time
    println!("{}\n{}", k, prover_time);

    let proof = transcript.finalize();

    // let strategy = SingleVerifier::new(&params_verifier);
    // let mut transcript = Blake2bRead::<_, _, Challenge255<_>>::init(&proof[..]);

    // verify_proof(
    //     &params_verifier,
    //     pk.get_vk(),
    //     strategy,
    //     &[&[]],
    //     &mut transcript,
    // )
    // .unwrap();
}

struct FuncCount {
    num_fft: usize, 
    num_extended_fft: usize,
    num_msm: usize,
    num_btree: usize,
}

fn dummy_proof<C: CurveAffine, ConcreteCircuit: Circuit<C::Scalar>>(
    params: &Params<C>,
    cs: &ConstraintSystem<C::Scalar>,
    domain: &EvaluationDomain<C::Scalar>,
) -> FuncCount {
    let l = 1;
    let mut num_fft = 0 as usize;
    let mut num_extended_fft = 0 as usize;
    let mut num_msm = 0 as usize;
    let mut num_btree = 0 as usize;

    // (instance, advice) calculate (poly, coset, commitment)

    // NOTE(sphere): ins_commit, pt += l * n_ins * commit_lagrange_t
    num_msm += l * cs.num_instance_columns;
    // NOTE(sphere): ins_poly, pt += l * n_ins + lagrange_to_coeff_t
    num_fft += l * cs.num_instance_columns;
    // NOTE(sphere): ins_coset, pt += l * n_ins + coeff_to_extended_t
    num_extended_fft += l * cs.num_instance_columns;
    // NOTE(sphere): adv_commit, pt += l * n_adv * commit_lagrange_t
    num_msm += l * cs.num_advice_columns;
    // NOTE(sphere): adv_poly, pt += l * n_adv * lagrange_to_coeff_t
    num_fft += l * cs.num_advice_columns;
    // NOTE(sphere): adv_coset, pt += l * n_adv * coeff_to_extended_t
    num_extended_fft += l * cs.num_advice_columns;

    // NOTE(sphere): pt += l * n_lookup * commit_permuted
    //      NOTE(sphere): BTree cost for A' and S'.
    let num_lookups = cs.lookups.len();
    num_btree += l * num_lookups;

    // Commit to permutations.
    // NOTE(sphere): l * perm_commit_t
    //      commit_lagrange: z
    let num_perm_slices = (cs.permutation.get_columns().len() + (cs.degree() - 1)) / (cs.degree() - 2);
    num_msm += num_perm_slices;
    //      lagrange_to_coeff: z
    num_fft += num_perm_slices;
    //      coeff_to_extended: z
    num_extended_fft += num_perm_slices;
    
    // NOTE(sphere): pt += lookup_commit_product
    //      commit_lagrange: z
    num_msm += num_lookups;
    //      lagrange_to_coeff: z
    num_fft += num_lookups;

    // Commit to the vanishing argument's random polynomial for blinding h(x_3)
    // NOTE(sphere): vanishing_commit
    //      commit: random_poly
    num_msm += 1;

    // Evaluate the h(X) polynomial
    // NOTE(sphere): evaluate_h 3 coeff_to_extended for each lookup argument
    num_extended_fft += 3 * num_lookups;

    // Construct the vanishing argument's h(X) commitments
    // NOTE(sphere): pt += vanishing_construct
    //      extended_to_coeff: h_poly
    num_extended_fft +=  1;
    //      commit: h_poly_i
    let num_h_pieces = ((domain.extended_len() as u64 + params.n - 1) / params.n) as usize;
    num_msm += num_h_pieces;

    // NOTE(sphere): evaluating ins / adv / fix only contains linear_terms.

    // NOTE(sphere): vanishing_evaluate only contains linear_terms.

    // NOTE(sphere): permutation_evaluate only contains linear_terms.

    // NOTE(sphere): permutation_construct_evaluate only contains linear_terms.

    // NOTE(sphere): lookups_evaluate only contains linear_terms.

    // NOTE(sphere): sum up number of queries.

    // NOTE(sphere): multiopen(shplonk).
    //      commit: h_x, h_x
    num_msm += 2;
    FuncCount {
        num_fft, 
        num_extended_fft,
        num_msm,
        num_btree,
    }
}


/// cost_model_main is to generate a main function to run the cost model for a circuit.
#[macro_export]
macro_rules! cost_model_main {
    ($cir:expr) => {
        use halo2_proofs::dev::{
            simulate_circuit,
            SimLinearResult,
            estimate,
        };

        fn main() {
            // NOTE(sphere): get k from args
            let mode = std::env::args().nth(1).expect("no running-mode given");
            let k = std::env::args().nth(2).expect("no circuit size given").parse().unwrap();
            // NOTE(sphere): estimate linear cost (cfg == simulate)
            let circuit = $cir;
            if mode.eq(&String::from("simulate")) {
                simulate_circuit::<Bn256, _>(circuit, k);
            } else if mode.eq(&String::from("estimate")) {
                // let res_path_1 = std::env::args().nth(3).expect("no circuit size given").parse().unwrap();
                // let res_path_2 = std::env::args().nth(4).expect("no circuit size given").parse().unwrap();
                // let res_1 = SimLinearResult::read(res_path_1);
                // let res_2 = SimLinearResult::read(res_path_2);
                let res_1 = SimLinearResult::new(10, 6292.0);
                let res_2 = SimLinearResult::new(14, 50092.0);
                let res = estimate::<Bn256, _>(circuit, res_1, res_2, k);
                res.print();
            } else {
                panic!("unrecognized format");
            }
        }
    }
}
