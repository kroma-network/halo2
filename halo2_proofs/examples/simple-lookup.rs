use halo2_proofs::arithmetic::FieldExt;
use halo2_proofs::circuit::{Layouter, SimpleFloorPlanner, Value};
use halo2_proofs::plonk::*;
use halo2_proofs::poly::commitment::ParamsProver;
use halo2_proofs::poly::kzg::commitment::{KZGCommitmentScheme, ParamsKZG};
use halo2_proofs::poly::kzg::multiopen::{ProverSHPLONK, VerifierSHPLONK};
use halo2_proofs::poly::kzg::strategy::SingleStrategy;
use halo2_proofs::poly::Rotation;
use halo2_proofs::transcript::{
    Blake2bRead, Blake2bWrite, Challenge255, TranscriptReadBuffer, TranscriptWriterBuffer,
};
use halo2curves::bn256::{Bn256, G1Affine};
use rand_core::SeedableRng;
use rand_xorshift::XorShiftRng;

use std::marker::PhantomData;

#[derive(Clone, Default)]
struct MyCircuit<F: FieldExt> {
    _marker: PhantomData<F>,
}

#[derive(Clone)]
struct MyConfig {
    selector: Selector,
    table: TableColumn,
    advice: Column<Advice>,
}

impl<F: FieldExt> Circuit<F> for MyCircuit<F> {
    type Config = MyConfig;
    type FloorPlanner = SimpleFloorPlanner;

    fn without_witnesses(&self) -> Self {
        Self::default()
    }

    fn configure(meta: &mut ConstraintSystem<F>) -> MyConfig {
        let config = MyConfig {
            selector: meta.complex_selector(),
            table: meta.lookup_table_column(),
            advice: meta.advice_column(),
        };

        meta.lookup("lookup", |meta| {
            let selector = meta.query_selector(config.selector);
            let not_selector = Expression::Constant(F::one()) - selector.clone();
            let advice = meta.query_advice(config.advice, Rotation::cur());
            vec![(selector * advice + not_selector, config.table)]
        });

        config
    }

    fn synthesize(&self, config: MyConfig, mut layouter: impl Layouter<F>) -> Result<(), Error> {
        layouter.assign_table(
            || "3-bit table",
            |mut table| {
                for row in 0u64..(1 << 3) {
                    table.assign_cell(
                        || format!("row {}", row),
                        config.table,
                        row as usize,
                        || Value::known(F::from(row + 1)),
                    )?;
                }

                Ok(())
            },
        )?;

        layouter.assign_region(
            || "assign values",
            |mut region| {
                for offset in 0u64..(1 << 4) {
                    config.selector.enable(&mut region, offset as usize)?;
                    region.assign_advice(
                        || format!("offset {}", offset),
                        config.advice,
                        offset as usize,
                        || Value::known(F::from((offset % 8) + 1)),
                    )?;
                }

                Ok(())
            },
        )
    }
}

fn main() {
    use halo2_proofs::dev::MockProver;
    use halo2curves::bn256::Fr;

    // ANCHOR: test-circuit
    // The number of rows in our circuit cannot exceed 2^k. Since our example
    // circuit is very small, we can pick a very small value here.
    let k = 5;

    // Prepare the private and public inputs to the circuit!

    // Instantiate the circuit with the private inputs.
    let circuit = MyCircuit::<Fr> {
        _marker: PhantomData,
    };

    // Given the correct public input, our circuit will verify.
    let s = Fr::from(2);
    let params = ParamsKZG::<Bn256>::unsafe_setup_with_s(k, s);
    let pk = keygen_pk2(&params, &circuit).expect("vk should not fail");
    // let vk = keygen_pk2(&params, &circuit).expect("vk should not fail");
    // let prover = MockProver::run(k, &circuit, vec![public_inputs.clone()]).unwrap();
    // assert_eq!(prover.verify(), Ok(()));

    // // If we try some other public input, the proof will fail!
    // public_inputs[0] += Fp::one();
    // let prover = MockProver::run(k, &circuit, vec![public_inputs]).unwrap();
    // assert!(prover.verify().is_err());
    // ANCHOR_END: test-circuit

    // let public_inputs = vec![];
    let public_inputs2: Vec<&[Fr]> = vec![];
    let public_inputs3 = vec![&public_inputs2[..]];

    let rng = XorShiftRng::from_seed([
        0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06, 0xbc,
        0xe5,
    ]);

    let proof = {
        let mut transcript = Blake2bWrite::<_, G1Affine, Challenge255<_>>::init(vec![]);

        create_proof::<KZGCommitmentScheme<Bn256>, ProverSHPLONK<_>, _, _, _, _>(
            &params,
            &pk,
            &[circuit],
            public_inputs3.as_slice(),
            rng,
            &mut transcript,
        )
        .expect("proof generation should not fail");

        transcript.finalize()
    };

    println!("{:?}", proof);

    let mut verifier_transcript = Blake2bRead::<_, G1Affine, Challenge255<_>>::init(&proof[..]);

    let verifier_params = params.verifier_params();
    let strategy = SingleStrategy::new(&params);
    verify_proof::<_, VerifierSHPLONK<_>, _, Blake2bRead<_, _, Challenge255<_>>, _>(
        &verifier_params,
        pk.get_vk(),
        strategy,
        public_inputs3.as_slice(),
        &mut verifier_transcript,
    )
    .unwrap();

    println!("Success");
}
