use std::marker::PhantomData;

use halo2_proofs::arithmetic::FieldExt;
use halo2_proofs::circuit::{Cell, Chip, Layouter, Region, SimpleFloorPlanner};
use halo2_proofs::dev::MockProver;
use halo2_proofs::plonk::*;
use halo2_proofs::poly::Rotation;
use pairing::bn256::Fr as Fp;

#[warn(dead_code)]
#[derive(Clone, Debug)]
struct MyConfig {
    lhs: Column<Advice>, 
    rhs: Column<Advice>, 
    out: Column<Advice>, 
    s_add: Selector, 
    q_even: Selector, 
    table_even: Column<Advice>,
    instance: Column<Instance>,
}

struct MyChip<F: FieldExt> {
    config: MyConfig,
    _marker: PhantomData<F>,
}

impl<F: FieldExt> Chip<F> for MyChip<F> {
    type Config = MyConfig;
    type Loaded = ();

    fn config(&self) -> &Self::Config {
        &self.config
    }

    fn loaded(&self) -> &Self::Loaded {
        &()
    }
}

/// The full circuit implementation.
///
/// In this struct we store the private input variables. We use `Option<F>` because
/// they won't have any value during key generation. During proving, if any of these
/// were `None` we would get an error.
#[derive(Default)]
struct MyCircuit<F: FieldExt> {
    out_lookup: Vec<F>, 
    x: Option<F>,
    d: Option<F>,
    t: usize, 
}

/// A variable representing a number.
#[derive(Clone)]
struct Number<F: FieldExt> {
    cell: Cell,
    value: Option<F>,
}

impl<F: FieldExt> MyChip<F> {
    fn construct(config: <Self as Chip<F>>::Config, _loaded: <Self as Chip<F>>::Loaded) -> Self {
        Self {
            config,
            _marker: PhantomData,
        }
    }

    fn load_private_lhs(
        &self,
        mut layouter: impl Layouter<F>,
        value: Option<F>,
    ) -> Result<Number<F>, Error> {
        let config = self.config();
        let mut num = None;
        layouter.assign_region(
            || "load private",
            |mut region| {
                let cell = region
                    .assign_advice(
                        || "private input",
                        config.lhs,
                        0,
                        || value.ok_or(Error::Synthesis),
                    )?
                    .cell();
                num = Some(Number { cell, value });
                Ok(())
            },
        )?;
        Ok(num.unwrap())
    }

    fn load_private_rhs(
        &self,
        mut layouter: impl Layouter<F>,
        value: Option<F>,
    ) -> Result<Number<F>, Error> {
        let config = self.config();
        let mut num = None;
        layouter.assign_region(
            || "load private",
            |mut region| {
                let cell = region
                    .assign_advice(
                        || "private input",
                        config.rhs,
                        0,
                        || value.ok_or(Error::Synthesis),
                    )?
                    .cell();
                num = Some(Number { cell, value });
                Ok(())
            },
        )?;
        Ok(num.unwrap())
    }

    fn add(
        &self,
        mut layouter: impl Layouter<F>,
        a: Number<F>,
        b: Number<F>,
    ) -> Result<Number<F>, Error> {
        let config = self.config();
        let mut out = None;
        layouter.assign_region(
            || "add",
            |mut region: Region<'_, F>| {

                config.s_add.enable(&mut region, 0)?;
                let lhs = region
                    .assign_advice(
                        || "lhs",
                        config.lhs,
                        0,
                        || a.value.ok_or(Error::Synthesis),
                    )?
                    .cell();
                let rhs = region
                    .assign_advice(
                        || "rhs",
                        config.rhs,
                        0,
                        || b.value.ok_or(Error::Synthesis),
                    )?
                    .cell();
                region.constrain_equal(a.cell, lhs)?;
                region.constrain_equal(b.cell, rhs)?;

                let value = a.value.and_then(|a| b.value.map(|b| a + b));
                let cell = region
                    .assign_advice(
                        || "lhs + rhs",
                        config.out,
                        0,
                        || value.ok_or(Error::Synthesis),
                    )?
                    .cell();

                config.q_even.enable(&mut region, 0)?;
                
                out = Some(Number { cell, value });
                Ok(())
            },
        )?;

        Ok(out.unwrap())
    }

    fn load_out_lookup(
        &self,
        mut layouter: impl Layouter<F>,
        values: &[F],
    ) -> Result<(), Error> {
        let config = self.config();
        layouter.assign_region(
            || "load values for even lookup table",
            |mut region| {
                for (offset, value) in values.iter().enumerate() {
                    region.assign_advice(
                        || "even table value",
                        config.table_even,
                        offset,
                        || Ok(*value),
                    )?;
                }

                Ok(())
            },
        )
    }

    fn expose_public(
        &self,
        mut layouter: impl Layouter<F>,
        num: Number<F>,
        row: usize,
    ) -> Result<(), Error> {
        layouter.constrain_instance(num.cell, self.config().instance, row)
    }
}

impl<F: FieldExt> Circuit<F> for MyCircuit<F> {
    type Config = MyConfig;
    type FloorPlanner = SimpleFloorPlanner;

    fn without_witnesses(&self) -> Self {
        Self::default()
    }

    fn configure(meta: &mut ConstraintSystem<F>) -> Self::Config {
        let lhs = meta.advice_column();
        let rhs = meta.advice_column();
        let out = meta.advice_column();
        let s_add = meta.selector();
        let q_even = meta.complex_selector();
        let table_even = meta.advice_column();
        let instance = meta.instance_column();

        meta.enable_equality(instance);
        meta.enable_equality(lhs);
        meta.enable_equality(rhs);
        meta.enable_equality(out);

        meta.create_gate("add", |meta| {
                let s_add = meta.query_selector(s_add);
                let lhs = meta.query_advice(lhs, Rotation::cur());
                let rhs = meta.query_advice(rhs, Rotation::cur());
                let out = meta.query_advice(out, Rotation::cur());
                vec![s_add * (lhs + rhs - out)]
        });

        meta.lookup_any("even number", |meta| {
            let q_even = meta.query_selector(q_even);
            let input = meta.query_advice(out, Rotation::cur());
            let table_even = meta.query_advice(table_even, Rotation::cur());

            vec![(q_even * input, table_even)]
        });

        MyConfig {
            lhs, 
            rhs, 
            out, 
            s_add, 
            q_even, 
            table_even,
            instance, 
        }
    }

    fn synthesize(
        &self,
        config: Self::Config,
        mut layouter: impl Layouter<F>,
    ) -> Result<(), Error> {
        let my_chip = MyChip::<F>::construct(config, ());
        my_chip.load_out_lookup(
            layouter.namespace(|| "witness even numbers"),
            &self.out_lookup, 
        )?;

        let mut x = my_chip.load_private_lhs(layouter.namespace(|| "load x"), self.x)?;
        let d = my_chip.load_private_rhs(layouter.namespace(|| "load d"), self.d)?;

        for _ in 0..self.t {
            let res = my_chip.add(layouter.namespace(|| "x + d"), x, d.clone())?;
            x = res;
        }
        my_chip.expose_public(layouter.namespace(|| "expose x"), x, 0)
    }
}

/// We build a circuit to calculate x + t * d. Suppose t = 2, Here is the layout. To 
/// include the lookup argument, we make a useless lookup argument on the out column. 
/// Since t = 2^k - 2, the total number of rows is 2^k - 1, which is the maximum under
/// the non-zero-knowledge setting.
/// | instance   | s_add | lhs   | rhs   | out       | q_even    | out_lookup    |
/// ------------------------------------------------------------------------------
/// | x + 2d     | 0     | x     | d     |           | 0         | x + d         |
/// |            | 1     | x     | d     | x + d     | 1         | x + 2d        |
/// |            | 1     | x + d | d     | x + 2d    | 1         | 0             |
/// ------------------------------------------------------------------------------
#[cfg(not(feature = "zero-knowledge"))]
fn main() {
    let k = 5;

    // Prepare the private and public inputs to the circuit!
    let x = Fp::zero();
    let d = Fp::from(2);

    // The maximum number of rows (2^k - 1) - 1 (load_private_input row).
    let t = (1 << k) - 2 as usize;  
    let y = x + d * Fp::from(t as u64);


    let mut out_lookup = Vec::with_capacity(t);
    for i in 0..t {
        out_lookup.push(x + d * Fp::from((i + 1) as u64));
    }

    // Instantiate the circuit with the private inputs.
    let circuit = MyCircuit {
        out_lookup, 
        x: Some(x),
        d: Some(d), 
        t, 
    };

    // Arrange the public input. We expose the result in row 0
    // of the instance column, so we position it there in our public inputs.
    let public_inputs = vec![y];

    // Given the correct public input, our circuit will verify.
    let prover = MockProver::run(k, &circuit, vec![public_inputs.clone()]).unwrap();
    assert_eq!(prover.verify(), Ok(()));

    // If we try with one more operation, the number of rows isn't enough!
    let t = (1 << k) - 1 as usize;  

    // Prepare the private and public inputs to the circuit!
    let x = Fp::zero();
    let d = Fp::from(2);

    let y = x + d * Fp::from(t as u64);

    let mut out_lookup = Vec::with_capacity(t);
    for i in 0..t {
        out_lookup.push(x + d * Fp::from((i + 1) as u64));
    }

    let circuit = MyCircuit {
        out_lookup, 
        x: Some(x),
        d: Some(d), 
        t, 
    };

    let public_inputs = vec![y];
    assert!(matches!(
        MockProver::<Fp>::run(k, &circuit, vec![public_inputs]).unwrap_err(),
        Error::NotEnoughRowsAvailable {
            current_k,
        } if current_k == k,
    ));
}

#[cfg(feature = "zero-knowledge")]
fn main() {}
