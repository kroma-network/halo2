use ff::BatchInvert;
use halo2_proofs::{
    arithmetic::FieldExt,
    circuit::{Layouter, SimpleFloorPlanner, Value},
    halo2curves::bn256::{Bn256, Fr},
    plonk::*,
    poly::{
        kzg::{
            commitment::{KZGCommitmentScheme, ParamsKZG},
            multiopen::{ProverSHPLONK, VerifierSHPLONK},
            strategy::SingleStrategy,
        },
        Rotation,
    },
    transcript::{
        Blake2bRead, Blake2bWrite, Challenge255, TranscriptReadBuffer, TranscriptWriterBuffer,
    },
};
use rand_core::{RngCore, SeedableRng};
use rand_xorshift::XorShiftRng;
use std::iter;

fn rand_2d_array<F: FieldExt, R: RngCore, const W: usize, const H: usize>(
    rng: &mut R,
) -> [[F; H]; W] {
    [(); W].map(|_| [(); H].map(|_| F::random(&mut *rng)))
}

fn shuffled<F: FieldExt, R: RngCore, const W: usize, const H: usize>(
    original: [[F; H]; W],
    rng: &mut R,
) -> [[F; H]; W] {
    let mut shuffled = original;

    for row in (1..H).rev() {
        let rand_row = (rng.next_u32() as usize) % row;
        for column in shuffled.iter_mut() {
            column.swap(row, rand_row);
        }
    }

    shuffled
}

#[derive(Clone)]
struct MyConfig<const W: usize> {
    q_shuffle: Selector,
    q_first: Selector,
    q_last: Selector,
    original: [Column<Advice>; W],
    shuffled: [Column<Advice>; W],
    theta: Challenge,
    gamma: Challenge,
    z: Column<Advice>,
}

impl<const W: usize> MyConfig<W> {
    fn configure<F: FieldExt>(meta: &mut ConstraintSystem<F>) -> Self {
        let [q_shuffle, q_first, q_last] = [(); 3].map(|_| meta.selector());
        // First phase
        let original = [(); W].map(|_| meta.advice_column_in(FirstPhase));
        let shuffled = [(); W].map(|_| meta.advice_column_in(FirstPhase));
        let [theta, gamma] = [(); 2].map(|_| meta.challenge_usable_after(FirstPhase));
        // Second phase
        let z = meta.advice_column_in(SecondPhase);

        meta.create_gate("z should start with 1", |meta| {
            let q_first = meta.query_selector(q_first);
            let z = meta.query_advice(z, Rotation::cur());
            let one = Expression::Constant(F::one());

            vec![q_first * (one - z)]
        });

        meta.create_gate("z should end with 1", |meta| {
            let q_last = meta.query_selector(q_last);
            let z = meta.query_advice(z, Rotation::cur());
            let one = Expression::Constant(F::one());

            vec![q_last * (one - z)]
        });

        meta.create_gate("z should have valid transition", |meta| {
            let q_shuffle = meta.query_selector(q_shuffle);
            let original = original.map(|advice| meta.query_advice(advice, Rotation::cur()));
            let shuffled = shuffled.map(|advice| meta.query_advice(advice, Rotation::cur()));
            let [theta, gamma] = [theta, gamma].map(|challenge| meta.query_challenge(challenge));
            let [z, z_w] =
                [Rotation::cur(), Rotation::next()].map(|rotation| meta.query_advice(z, rotation));

            // Compress
            let original = original
                .iter()
                .cloned()
                .reduce(|acc, a| acc * theta.clone() + a)
                .unwrap();
            let shuffled = shuffled
                .iter()
                .cloned()
                .reduce(|acc, a| acc * theta.clone() + a)
                .unwrap();

            vec![q_shuffle * (z * (original + gamma.clone()) - z_w * (shuffled + gamma))]
        });

        Self {
            q_shuffle,
            q_first,
            q_last,
            original,
            shuffled,
            theta,
            gamma,
            z,
        }
    }
}

#[derive(Clone, Default)]
struct MyCircuit<F: FieldExt, const W: usize, const H: usize> {
    original: Value<[[F; H]; W]>,
    shuffled: Value<[[F; H]; W]>,
}

impl<F: FieldExt, const W: usize, const H: usize> MyCircuit<F, W, H> {
    fn rand<R: RngCore>(rng: &mut R) -> Self {
        let original = rand_2d_array::<F, _, W, H>(rng);
        let shuffled = shuffled(original, rng);

        Self {
            original: Value::known(original),
            shuffled: Value::known(shuffled),
        }
    }
}

impl<F: FieldExt, const W: usize, const H: usize> Circuit<F> for MyCircuit<F, W, H> {
    type Config = MyConfig<W>;
    type FloorPlanner = SimpleFloorPlanner;

    fn without_witnesses(&self) -> Self {
        Self::default()
    }

    fn configure(meta: &mut ConstraintSystem<F>) -> Self::Config {
        MyConfig::configure(meta)
    }

    fn synthesize(
        &self,
        config: Self::Config,
        mut layouter: impl Layouter<F>,
    ) -> Result<(), Error> {
        let theta = layouter.get_challenge(config.theta);
        let gamma = layouter.get_challenge(config.gamma);

        layouter.assign_region(
            || "Shuffle original into shuffled",
            |mut region| {
                // Keygen
                config.q_first.enable(&mut region, 0)?;
                config.q_last.enable(&mut region, H)?;
                for offset in 0..H {
                    config.q_shuffle.enable(&mut region, offset)?;
                }

                // First phase
                for (idx, (&column, values)) in config
                    .original
                    .iter()
                    .zip(self.original.transpose_array().iter())
                    .enumerate()
                {
                    for (offset, &value) in values.transpose_array().iter().enumerate() {
                        region.assign_advice(
                            || format!("original[{}][{}]", idx, offset),
                            column,
                            offset,
                            || value,
                        )?;
                    }
                }
                for (idx, (&column, values)) in config
                    .shuffled
                    .iter()
                    .zip(self.shuffled.transpose_array().iter())
                    .enumerate()
                {
                    for (offset, &value) in values.transpose_array().iter().enumerate() {
                        region.assign_advice(
                            || format!("shuffled[{}][{}]", idx, offset),
                            column,
                            offset,
                            || value,
                        )?;
                    }
                }

                // Second phase
                let z = self.original.zip(self.shuffled).zip(theta).zip(gamma).map(
                    |(((original, shuffled), theta), gamma)| {
                        let mut product = vec![F::zero(); H];
                        for (idx, product) in product.iter_mut().enumerate() {
                            let mut compressed = F::zero();
                            for value in shuffled.iter() {
                                compressed *= theta;
                                compressed += value[idx];
                            }

                            *product = compressed + gamma
                        }

                        product.iter_mut().batch_invert();

                        for (idx, product) in product.iter_mut().enumerate() {
                            let mut compressed = F::zero();
                            for value in original.iter() {
                                compressed *= theta;
                                compressed += value[idx];
                            }

                            *product *= compressed + gamma
                        }

                        #[allow(clippy::let_and_return)]
                        let z = iter::once(F::one())
                            .chain(product)
                            .scan(F::one(), |state, cur| {
                                *state *= &cur;
                                Some(*state)
                            })
                            .collect::<Vec<_>>();

                        #[cfg(feature = "sanity-checks")]
                        assert_eq!(F::one(), *z.last().unwrap());

                        z
                    },
                );
                for (offset, value) in z.transpose_vec(H + 1).into_iter().enumerate() {
                    region.assign_advice(
                        || format!("z[{}]", offset),
                        config.z,
                        offset,
                        || value,
                    )?;
                }

                Ok(())
            },
        )
    }
}

fn main() {
    const W: usize = 2;
    const H: usize = 8;
    const K: u32 = 4;

    let mut rng_for_table = XorShiftRng::from_seed([
        0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06, 0xbc,
        0xe5,
    ]);
    let circuit = MyCircuit::<Fr, W, H>::rand(&mut rng_for_table);

    let s = Fr::from(2);
    let params = ParamsKZG::<Bn256>::unsafe_setup_with_s(K, s);
    let pk = keygen_pk2(&params, &circuit).unwrap();

    let rng = XorShiftRng::from_seed([
        0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06, 0xbc,
        0xe5,
    ]);

    let proof = {
        let mut transcript = Blake2bWrite::<_, _, Challenge255<_>>::init(vec![]);

        create_proof::<KZGCommitmentScheme<Bn256>, ProverSHPLONK<_>, _, _, _, _>(
            &params,
            &pk,
            &[circuit],
            &[&[]],
            rng,
            &mut transcript,
        )
        .expect("proof generation should not fail");

        transcript.finalize()
    };
    println!("proof: \n{:?}", proof);

    let strategy = SingleStrategy::new(&params);
    let mut transcript = Blake2bRead::<_, _, Challenge255<_>>::init(&proof[..]);
    verify_proof::<KZGCommitmentScheme<Bn256>, VerifierSHPLONK<Bn256>, _, _, _>(
        &params,
        pk.get_vk(),
        strategy,
        &[&[], &[]],
        &mut transcript,
    )
    .unwrap();
}
