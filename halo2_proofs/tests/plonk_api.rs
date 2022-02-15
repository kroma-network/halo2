#![allow(clippy::many_single_char_names)]
#![allow(clippy::op_ref)]

use assert_matches::assert_matches;
use halo2_proofs::arithmetic::{CurveAffine, FieldExt};
use halo2_proofs::circuit::{Cell, Layouter, SimpleFloorPlanner, Value};
use halo2_proofs::dev::MockProver;
use halo2_proofs::plonk::{
    create_proof, keygen_pk, keygen_vk, verify_proof, Advice, Assigned, BatchVerifier, Circuit,
    Column, ConstraintSystem, Error, Fixed, SingleVerifier, TableColumn, VerifyingKey,
};

use halo2_proofs::poly::{
    commitment::{Params, ParamsVerifier},
    Rotation,
};
use halo2_proofs::transcript::{Blake2bRead, Blake2bWrite, Challenge255};
use rand_core::OsRng;
use std::marker::PhantomData;

use pairing::bn256::Fr as Fp;
use pairing::bn256::{Bn256, G1Affine};

#[test]
fn plonk_api() {
    const K: u32 = 5;
    let public_inputs_size = 1;

    /// This represents an advice column at a certain row in the ConstraintSystem
    #[derive(Copy, Clone, Debug)]
    pub struct Variable(Column<Advice>, usize);

    // Initialize the polynomial commitment parameters
    let params: Params<G1Affine> = Params::<G1Affine>::unsafe_setup::<Bn256>(K);
    let params_verifier: ParamsVerifier<Bn256> = params.verifier(public_inputs_size).unwrap();

    #[derive(Clone)]
    struct PlonkConfig {
        a: Column<Advice>,
        b: Column<Advice>,
        c: Column<Advice>,
        d: Column<Advice>,
        e: Column<Advice>,

        sa: Column<Fixed>,
        sb: Column<Fixed>,
        sc: Column<Fixed>,
        sm: Column<Fixed>,
        sp: Column<Fixed>,
        sl: TableColumn,
    }

    #[allow(clippy::type_complexity)]
    trait StandardCs<FF: FieldExt> {
        fn raw_multiply<F>(
            &self,
            layouter: &mut impl Layouter<FF>,
            f: F,
        ) -> Result<(Cell, Cell, Cell), Error>
        where
            F: FnMut() -> Value<(Assigned<FF>, Assigned<FF>, Assigned<FF>)>;
        fn raw_add<F>(
            &self,
            layouter: &mut impl Layouter<FF>,
            f: F,
        ) -> Result<(Cell, Cell, Cell), Error>
        where
            F: FnMut() -> Value<(Assigned<FF>, Assigned<FF>, Assigned<FF>)>;
        fn copy(&self, layouter: &mut impl Layouter<FF>, a: Cell, b: Cell) -> Result<(), Error>;
        fn public_input<F>(&self, layouter: &mut impl Layouter<FF>, f: F) -> Result<Cell, Error>
        where
            F: FnMut() -> Value<FF>;
        fn lookup_table(
            &self,
            layouter: &mut impl Layouter<FF>,
            values: &[FF],
        ) -> Result<(), Error>;
    }

    #[derive(Clone)]
    struct MyCircuit<F: FieldExt> {
        a: Value<F>,
        lookup_table: Vec<F>,
    }

    struct StandardPlonk<F: FieldExt> {
        config: PlonkConfig,
        _marker: PhantomData<F>,
    }

    impl<FF: FieldExt> StandardPlonk<FF> {
        fn new(config: PlonkConfig) -> Self {
            StandardPlonk {
                config,
                _marker: PhantomData,
            }
        }
    }

    impl<FF: FieldExt> StandardCs<FF> for StandardPlonk<FF> {
        fn raw_multiply<F>(
            &self,
            layouter: &mut impl Layouter<FF>,
            mut f: F,
        ) -> Result<(Cell, Cell, Cell), Error>
        where
            F: FnMut() -> Value<(Assigned<FF>, Assigned<FF>, Assigned<FF>)>,
        {
            layouter.assign_region(
                || "raw_multiply",
                |mut region| {
                    let mut value = None;
                    let lhs = region.assign_advice(
                        || "lhs",
                        self.config.a,
                        0,
                        || {
                            value = Some(f());
                            value.unwrap().map(|v| v.0)
                        },
                    )?;
                    region.assign_advice(
                        || "lhs^4",
                        self.config.d,
                        0,
                        || value.unwrap().map(|v| v.0).square().square(),
                    )?;
                    let rhs = region.assign_advice(
                        || "rhs",
                        self.config.b,
                        0,
                        || value.unwrap().map(|v| v.1),
                    )?;
                    region.assign_advice(
                        || "rhs^4",
                        self.config.e,
                        0,
                        || value.unwrap().map(|v| v.1).square().square(),
                    )?;
                    let out = region.assign_advice(
                        || "out",
                        self.config.c,
                        0,
                        || value.unwrap().map(|v| v.2),
                    )?;

                    region.assign_fixed(|| "a", self.config.sa, 0, || Value::known(FF::zero()))?;
                    region.assign_fixed(|| "b", self.config.sb, 0, || Value::known(FF::zero()))?;
                    region.assign_fixed(|| "c", self.config.sc, 0, || Value::known(FF::one()))?;
                    region.assign_fixed(
                        || "a * b",
                        self.config.sm,
                        0,
                        || Value::known(FF::one()),
                    )?;
                    Ok((lhs.cell(), rhs.cell(), out.cell()))
                },
            )
        }
        fn raw_add<F>(
            &self,
            layouter: &mut impl Layouter<FF>,
            mut f: F,
        ) -> Result<(Cell, Cell, Cell), Error>
        where
            F: FnMut() -> Value<(Assigned<FF>, Assigned<FF>, Assigned<FF>)>,
        {
            layouter.assign_region(
                || "raw_add",
                |mut region| {
                    let mut value = None;
                    let lhs = region.assign_advice(
                        || "lhs",
                        self.config.a,
                        0,
                        || {
                            value = Some(f());
                            value.unwrap().map(|v| v.0)
                        },
                    )?;
                    region.assign_advice(
                        || "lhs^4",
                        self.config.d,
                        0,
                        || value.unwrap().map(|v| v.0).square().square(),
                    )?;
                    let rhs = region.assign_advice(
                        || "rhs",
                        self.config.b,
                        0,
                        || value.unwrap().map(|v| v.1),
                    )?;
                    region.assign_advice(
                        || "rhs^4",
                        self.config.e,
                        0,
                        || value.unwrap().map(|v| v.1).square().square(),
                    )?;
                    let out = region.assign_advice(
                        || "out",
                        self.config.c,
                        0,
                        || value.unwrap().map(|v| v.2),
                    )?;

                    region.assign_fixed(|| "a", self.config.sa, 0, || Value::known(FF::one()))?;
                    region.assign_fixed(|| "b", self.config.sb, 0, || Value::known(FF::one()))?;
                    region.assign_fixed(|| "c", self.config.sc, 0, || Value::known(FF::one()))?;
                    region.assign_fixed(
                        || "a * b",
                        self.config.sm,
                        0,
                        || Value::known(FF::zero()),
                    )?;
                    Ok((lhs.cell(), rhs.cell(), out.cell()))
                },
            )
        }
        fn copy(
            &self,
            layouter: &mut impl Layouter<FF>,
            left: Cell,
            right: Cell,
        ) -> Result<(), Error> {
            layouter.assign_region(
                || "copy",
                |mut region| {
                    region.constrain_equal(left, right)?;
                    region.constrain_equal(left, right)
                },
            )
        }
        fn public_input<F>(&self, layouter: &mut impl Layouter<FF>, mut f: F) -> Result<Cell, Error>
        where
            F: FnMut() -> Value<FF>,
        {
            layouter.assign_region(
                || "public_input",
                |mut region| {
                    let value = region.assign_advice(|| "value", self.config.a, 0, &mut f)?;
                    region.assign_fixed(
                        || "public",
                        self.config.sp,
                        0,
                        || Value::known(FF::one()),
                    )?;

                    Ok(value.cell())
                },
            )
        }
        fn lookup_table(
            &self,
            layouter: &mut impl Layouter<FF>,
            values: &[FF],
        ) -> Result<(), Error> {
            layouter.assign_table(
                || "",
                |mut table| {
                    for (index, &value) in values.iter().enumerate() {
                        table.assign_cell(
                            || "table col",
                            self.config.sl,
                            index,
                            || Value::known(value),
                        )?;
                    }
                    Ok(())
                },
            )?;
            Ok(())
        }
    }

    impl<F: FieldExt> Circuit<F> for MyCircuit<F> {
        type Config = PlonkConfig;
        type FloorPlanner = SimpleFloorPlanner;

        fn without_witnesses(&self) -> Self {
            Self {
                a: Value::unknown(),
                lookup_table: self.lookup_table.clone(),
            }
        }

        fn configure(meta: &mut ConstraintSystem<F>) -> PlonkConfig {
            let e = meta.advice_column();
            let a = meta.advice_column();
            let b = meta.advice_column();
            let sf = meta.fixed_column();
            let c = meta.advice_column();
            let d = meta.advice_column();
            let p = meta.instance_column();

            meta.enable_equality(a);
            meta.enable_equality(b);
            meta.enable_equality(c);

            let sm = meta.fixed_column();
            let sa = meta.fixed_column();
            let sb = meta.fixed_column();
            let sc = meta.fixed_column();
            let sp = meta.fixed_column();
            let sl = meta.lookup_table_column();

            /*
             *   A         B      ...  sl
             * [
             *   instance  0      ...  0
             *   a         a      ...  0
             *   a         a^2    ...  0
             *   a         a      ...  0
             *   a         a^2    ...  0
             *   ...       ...    ...  ...
             *   ...       ...    ...  instance
             *   ...       ...    ...  a
             *   ...       ...    ...  a
             *   ...       ...    ...  0
             * ]
             */

            meta.lookup("lookup", |meta| {
                let a_ = meta.query_any(a, Rotation::cur());
                vec![(a_, sl)]
            });

            meta.create_gate("Combined add-mult", |meta| {
                let d = meta.query_advice(d, Rotation::next());
                let a = meta.query_advice(a, Rotation::cur());
                let sf = meta.query_fixed(sf, Rotation::cur());
                let e = meta.query_advice(e, Rotation::prev());
                let b = meta.query_advice(b, Rotation::cur());
                let c = meta.query_advice(c, Rotation::cur());

                let sa = meta.query_fixed(sa, Rotation::cur());
                let sb = meta.query_fixed(sb, Rotation::cur());
                let sc = meta.query_fixed(sc, Rotation::cur());
                let sm = meta.query_fixed(sm, Rotation::cur());

                vec![a.clone() * sa + b.clone() * sb + a * b * sm - (c * sc) + sf * (d * e)]
            });

            meta.create_gate("Public input", |meta| {
                let a = meta.query_advice(a, Rotation::cur());
                let p = meta.query_instance(p, Rotation::cur());
                let sp = meta.query_fixed(sp, Rotation::cur());

                vec![sp * (a - p)]
            });

            meta.enable_equality(sf);
            meta.enable_equality(e);
            meta.enable_equality(d);
            meta.enable_equality(p);
            meta.enable_equality(sm);
            meta.enable_equality(sa);
            meta.enable_equality(sb);
            meta.enable_equality(sc);
            meta.enable_equality(sp);

            PlonkConfig {
                a,
                b,
                c,
                d,
                e,
                sa,
                sb,
                sc,
                sm,
                sp,
                sl,
            }
        }

        fn synthesize(
            &self,
            config: PlonkConfig,
            mut layouter: impl Layouter<F>,
        ) -> Result<(), Error> {
            let cs = StandardPlonk::new(config);

            let _ = cs.public_input(&mut layouter, || Value::known(F::one() + F::one()))?;

            for _ in 0..10 {
                let a: Value<Assigned<_>> = self.a.into();
                let mut a_squared = Value::unknown();
                let (a0, _, c0) = cs.raw_multiply(&mut layouter, || {
                    a_squared = a.square();
                    a.zip(a_squared).map(|(a, a_squared)| (a, a, a_squared))
                })?;
                let (a1, b1, _) = cs.raw_add(&mut layouter, || {
                    let fin = a_squared + a;
                    a.zip(a_squared)
                        .zip(fin)
                        .map(|((a, a_squared), fin)| (a, a_squared, fin))
                })?;
                cs.copy(&mut layouter, a0, a1)?;
                cs.copy(&mut layouter, b1, c0)?;
            }

            cs.lookup_table(&mut layouter, &self.lookup_table)?;

            Ok(())
        }
    }

    let a = Fp::from(2834758237) * Fp::ZETA;
    let instance = Fp::one() + Fp::one();
    let lookup_table = vec![instance, a, a, Fp::zero()];

    let empty_circuit: MyCircuit<Fp> = MyCircuit {
        a: Value::unknown(),
        lookup_table: lookup_table.clone(),
    };

    let circuit: MyCircuit<Fp> = MyCircuit {
        a: Value::known(a),
        lookup_table,
    };

    // Check that we get an error if we try to initialize the proving key with a value of
    // k that is too small for the minimum required number of rows.
    let much_too_small_params: Params<G1Affine> = Params::<G1Affine>::unsafe_setup::<Bn256>(1);
    assert_matches!(
        keygen_vk(&much_too_small_params, &empty_circuit),
        Err(Error::NotEnoughRowsAvailable {
            current_k,
        }) if current_k == 1
    );

    // Check that we get an error if we try to initialize the proving key with a value of
    // k that is too small for the number of rows the circuit uses.
    let slightly_too_small_params: Params<G1Affine> =
        Params::<G1Affine>::unsafe_setup::<Bn256>(K - 1);
    assert_matches!(
        keygen_vk(&slightly_too_small_params, &empty_circuit),
        Err(Error::NotEnoughRowsAvailable {
            current_k,
        }) if current_k == K - 1
    );

    // Initialize the proving key
    let vk = keygen_vk(&params, &empty_circuit).expect("keygen_vk should not fail");
    let pk = keygen_pk(&params, vk, &empty_circuit).expect("keygen_pk should not fail");

    let pubinputs = vec![instance];

    // Check this circuit is satisfied.
    let prover = match MockProver::run(K, &circuit, vec![pubinputs.clone()]) {
        Ok(prover) => prover,
        Err(e) => panic!("{:?}", e),
    };
    assert_eq!(prover.verify(), Ok(()));

    if std::env::var_os("HALO2_PLONK_TEST_GENERATE_NEW_PROOF").is_some() {
        let mut transcript = Blake2bWrite::<_, _, Challenge255<_>>::init(vec![]);
        // Create a proof
        create_proof(
            &params,
            &pk,
            &[circuit.clone(), circuit.clone()],
            &[&[&[instance]], &[&[instance]]],
            OsRng,
            &mut transcript,
        )
        .expect("proof generation should not fail");
        let proof: Vec<u8> = transcript.finalize();

        std::fs::write("plonk_api_proof.bin", &proof[..])
            .expect("should succeed to write new proof");
    }

    {
        // Check that a hardcoded proof is satisfied
        let proof = include_bytes!("plonk_api_proof.bin");
        let strategy = SingleVerifier::new(&params);
        let mut transcript = Blake2bRead::<_, _, Challenge255<_>>::init(&proof[..]);
        assert!(verify_proof(
            &params,
            pk.get_vk(),
            strategy,
            &[&[&pubinputs[..]], &[&pubinputs[..]]],
            &mut transcript,
        )
        .is_ok());
    }

    for _ in 0..10 {
        let mut transcript = Blake2bWrite::<_, _, Challenge255<_>>::init(vec![]);
        // Create a proof
        create_proof(
            &params,
            &pk,
            &[circuit.clone(), circuit.clone()],
            &[&[&[instance]], &[&[instance]]],
            OsRng,
            &mut transcript,
        )
        .expect("proof generation should not fail");
        let proof: Vec<u8> = transcript.finalize();

        // Test single-verifier strategy.
        {
            let strategy = SingleVerifier::new(&params_verifier);
            let mut transcript = Blake2bRead::<_, _, Challenge255<_>>::init(&proof[..]);
            assert!(verify_proof(
                &params_verifier,
                pk.get_vk(),
                strategy,
                &[&[&pubinputs[..]], &[&pubinputs[..]]],
                &mut transcript,
            )
            .is_ok());
        }

        //
        // Test batch-verifier strategy.
        //

        {
            let strategy = BatchVerifier::new(&params_verifier, OsRng);

            // First proof.
            let mut transcript = Blake2bRead::<_, _, Challenge255<_>>::init(&proof[..]);
            let strategy = verify_proof(
                &params_verifier,
                pk.get_vk(),
                strategy,
                &[&[&pubinputs[..]], &[&pubinputs[..]]],
                &mut transcript,
            )
            .unwrap();

            // Write and then read the verification key in between (to check round-trip
            // serialization).
            // TODO: Figure out whether https://github.com/zcash/halo2/issues/449 should
            // be caught by this, or if it is caused by downstream changes to halo2.
            let mut vk_buffer = vec![];
            pk.get_vk().write(&mut vk_buffer).unwrap();
            let vk =
                VerifyingKey::<G1Affine>::read::<_, MyCircuit<Fp>>(&mut &vk_buffer[..], &params)
                    .unwrap();

            // "Second" proof (just the first proof again).
            let mut transcript = Blake2bRead::<_, _, Challenge255<_>>::init(&proof[..]);
            let strategy = verify_proof(
                &params_verifier,
                &vk,
                strategy,
                &[&[&pubinputs[..]], &[&pubinputs[..]]],
                &mut transcript,
            )
            .unwrap();

            // Check the batch.
            assert!(batch.finalize(&params, pk.get_vk()));
        }
    }
}
