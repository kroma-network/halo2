use std::marker::PhantomData;

use ff::{Field, PrimeField};
use halo2_proofs::{
    bn254::Blake2bWrite,
    circuit::{Layouter, SimpleFloorPlanner, Value},
    consts::SEED,
    plonk::{
        create_proof, keygen_pk2, Advice, Circuit, Column, ConstraintSystem, Error, Expression,
        Selector, TableColumn,
    },
    poly::{
        kzg::{
            commitment::{KZGCommitmentScheme, ParamsKZG},
            multiopen::ProverGWC,
        },
        Rotation,
    },
    transcript::{Challenge255, TranscriptWriterBuffer},
    xor_shift_rng::XORShiftRng,
};
use halo2curves::bn256::{Bn256, G1Affine};
use rand_core::SeedableRng;

#[derive(Clone, Default)]
struct SimpleLookupCircuit<F: Field> {
    _marker: PhantomData<F>,
}

#[derive(Clone)]
struct SimpleLookupConfig {
    selector: Selector,
    table: TableColumn,
    advice: Column<Advice>,
}

impl<F: PrimeField> Circuit<F> for SimpleLookupCircuit<F> {
    type Config = SimpleLookupConfig;
    type FloorPlanner = SimpleFloorPlanner;

    fn without_witnesses(&self) -> Self {
        Self::default()
    }

    fn configure(meta: &mut ConstraintSystem<F>) -> SimpleLookupConfig {
        let config = SimpleLookupConfig {
            selector: meta.complex_selector(),
            table: meta.lookup_table_column(),
            advice: meta.advice_column(),
        };

        meta.lookup("lookup", |meta| {
            let selector = meta.query_selector(config.selector);
            let not_selector = Expression::Constant(F::ONE) - selector.clone();
            let advice = meta.query_advice(config.advice, Rotation::cur());
            vec![(selector * advice + not_selector, config.table)]
        });

        config
    }

    fn synthesize(
        &self,
        config: SimpleLookupConfig,
        mut layouter: impl Layouter<F>,
    ) -> Result<(), Error> {
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
    let vec = vec![
        1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1,
    ];
    let a = vec.binary_search(&1);
    println!("{:?}", a);

    use halo2curves::bn256::Fr;

    env_logger::init();

    // ANCHOR: test-circuit
    // The number of rows in our circuit cannot exceed 2^k. Since our example
    // circuit is very small, we can pick a very small value here.
    let k = 5;

    // Instantiate the circuit with the private inputs.
    let circuit = SimpleLookupCircuit::<Fr> {
        _marker: PhantomData,
    };
    // Arrange the public input.
    let public_inputs = vec![];
    let public_inputs2 = vec![&public_inputs[..], &public_inputs[..]];

    let s = Fr::from(2);
    let params = ParamsKZG::<Bn256>::unsafe_setup_with_s(k, s.clone());
    let pk = keygen_pk2(&params, &circuit).expect("vk should not fail");

    let rng = XORShiftRng::from_seed(SEED);

    let mut transcript = Blake2bWrite::<_, G1Affine, Challenge255<_>>::init(vec![]);

    create_proof::<KZGCommitmentScheme<Bn256>, ProverGWC<_>, _, _, _, _>(
        &params,
        &pk,
        &[circuit.clone(), circuit.clone()],
        public_inputs2.as_slice(),
        rng.clone(),
        &mut transcript,
    )
    .expect("proof generation should not fail");

    let proof = transcript.finalize();

    println!("done!");
}
