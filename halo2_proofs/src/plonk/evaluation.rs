use crate::multicore;
use crate::plonk::lookup::prover::Committed;
use crate::plonk::permutation::Argument;
use crate::plonk::{circuit, lookup, permutation, Any, ProvingKey};
use crate::poly::Basis;
use crate::{
    arithmetic::{eval_polynomial, parallelize, CurveAffine, FieldExt},
    poly::{
        commitment::Params, Coeff, EvaluationDomain, ExtendedLagrangeCoeff, LagrangeCoeff,
        Polynomial, ProverQuery, Rotation,
    },
    transcript::{EncodedChallenge, TranscriptWrite},
};
use circuit::Column;
use group::prime::PrimeCurve;
use group::{
    ff::{BatchInvert, Field},
    Curve,
};
use itertools::Itertools;
use std::any::TypeId;
use std::collections::BTreeSet;
use std::convert::TryInto;
use std::num::ParseIntError;
use std::ops::IndexMut;
use std::slice;
use std::{
    collections::BTreeMap,
    iter,
    ops::{Index, Mul, MulAssign},
};

use super::{ConstraintSystem, Expression};

/// Return the index in the polynomial of size `isize` after rotation `rot`.
fn get_rotation_idx(idx: usize, rot: i32, rot_scale: i32, isize: i32) -> usize {
    (((idx as i32) + (rot * rot_scale)).rem_euclid(isize)) as usize
}

/// Value used in a calculation
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd)]
pub enum ValueSource {
    /// This is a constant value
    Constant(usize),
    /// This is an intermediate value
    Intermediate(usize),
    /// This is a fixed column
    Fixed(usize, usize),
    /// This is an advice (witness) column
    Advice(usize, usize),
    /// This is an instance (external) column
    Instance(usize, usize),
    /// beta
    Beta(),
    /// gamma
    Gamma(),
    /// theta
    Theta(),
    /// y
    Y(),
    /// Previous value
    PreviousValue(),
}

impl Default for ValueSource {
    fn default() -> Self {
        ValueSource::Constant(0)
    }
}

impl ValueSource {
    /// Get the value for this source
    pub fn get<F: Field, B: Basis>(
        &self,
        rotations: &[usize],
        constants: &[F],
        intermediates: &[F],
        fixed_values: &[Polynomial<F, B>],
        advice_values: &[Option<Polynomial<F, B>>],
        instance_values: &[Option<Polynomial<F, B>>],
        beta: &F,
        gamma: &F,
        theta: &F,
        y: &F,
        previous_value: &F,
    ) -> F {
        match self {
            ValueSource::Constant(idx) => constants[*idx],
            ValueSource::Intermediate(idx) => intermediates[*idx],
            ValueSource::Fixed(column_index, rotation) => {
                fixed_values[*column_index][rotations[*rotation]]
            }
            ValueSource::Advice(column_index, rotation) => {
                advice_values[*column_index].as_ref().unwrap()[rotations[*rotation]]
            }
            ValueSource::Instance(column_index, rotation) => {
                instance_values[*column_index].as_ref().unwrap()[rotations[*rotation]]
            }
            ValueSource::Beta() => *beta,
            ValueSource::Gamma() => *gamma,
            ValueSource::Theta() => *theta,
            ValueSource::Y() => *y,
            ValueSource::PreviousValue() => *previous_value,
        }
    }
}

/// Calculation
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Calculation {
    /// This is an addition
    Add(ValueSource, ValueSource),
    /// This is a subtraction
    Sub(ValueSource, ValueSource),
    /// This is a product
    Mul(ValueSource, ValueSource),
    /// This is a square
    Square(ValueSource),
    /// This is a double
    Double(ValueSource),
    /// This is a negation
    Negate(ValueSource),
    /// This is Horner's rule: `val = a; val = val * c + b[]`
    Horner(ValueSource, Vec<ValueSource>, ValueSource),
    /// This is a simple assignment
    Store(ValueSource),
}

impl Calculation {
    /// Get the resulting value of this calculation
    pub fn evaluate<F: Field, B: Basis>(
        &self,
        rotations: &[usize],
        constants: &[F],
        intermediates: &[F],
        fixed_values: &[Polynomial<F, B>],
        advice_values: &[Option<Polynomial<F, B>>],
        instance_values: &[Option<Polynomial<F, B>>],
        beta: &F,
        gamma: &F,
        theta: &F,
        y: &F,
        previous_value: &F,
    ) -> F {
        let get_value = |value: &ValueSource| {
            value.get(
                rotations,
                constants,
                intermediates,
                fixed_values,
                advice_values,
                instance_values,
                beta,
                gamma,
                theta,
                y,
                previous_value,
            )
        };
        match self {
            Calculation::Add(a, b) => get_value(a) + get_value(b),
            Calculation::Sub(a, b) => get_value(a) - get_value(b),
            Calculation::Mul(a, b) => get_value(a) * get_value(b),
            Calculation::Square(v) => get_value(v).square(),
            Calculation::Double(v) => get_value(v).double(),
            Calculation::Negate(v) => -get_value(v),
            Calculation::Horner(start_value, parts, factor) => {
                let factor = get_value(factor);
                let mut value = get_value(start_value);
                for part in parts.iter() {
                    value = value * factor + get_value(part);
                }
                value
            }
            Calculation::Store(v) => get_value(v),
        }
    }
}

/// Evaluator
#[derive(Clone, Default, Debug)]
pub struct Evaluator<C: CurveAffine> {
    /// Used to control the alloc and release of coset polys
    pub mem_plan: MemPlan,
    ///  Custom gates evalution
    pub custom_gates: Vec<GraphEvaluator<C>>,
    ///  Lookups evalution with evaluation group id
    pub lookups: Vec<Vec<(GraphEvaluator<C>, usize)>>,
}

/// GraphEvaluator
#[derive(Clone, Debug)]
pub struct GraphEvaluator<C: CurveAffine> {
    /// Constants
    pub constants: Vec<C::ScalarExt>,
    /// Rotations
    pub rotations: Vec<i32>,
    /// Calculations
    pub calculations: Vec<CalculationInfo>,
    /// Number of intermediates
    pub num_intermediates: usize,
}

/// EvaluationData
#[derive(Default, Debug)]
pub struct EvaluationData<C: CurveAffine> {
    /// Intermediates
    pub intermediates: Vec<C::ScalarExt>,
    /// Rotations
    pub rotations: Vec<usize>,
}

/// CaluclationInfo
#[derive(Clone, Debug)]
pub struct CalculationInfo {
    /// Calculation
    pub calculation: Calculation,
    /// Target
    pub target: usize,
}

#[derive(Debug, Default, Clone)]
pub struct MemoryPlanOfEvalutionGroup {
    pub alloc_before: Vec<circuit::Column<Any>>,
    pub release_after: Vec<circuit::Column<Any>>,
}

#[derive(Default, Debug, Clone)]
pub struct MemPlan {
    pub gates: Vec<MemoryPlanOfEvalutionGroup>,
    pub lookups: Vec<MemoryPlanOfEvalutionGroup>,
}

impl MemPlan {
    // returns (fft_num, max_coset_poly_num)
    fn simulate_execution(&self) -> (usize, usize) {
        // simulate the execution of this coset memory plan
        let simulate = |plans: &[MemoryPlanOfEvalutionGroup]| {
            let mut coset_polys: BTreeSet<circuit::Column<Any>> = Default::default();
            let mut coset_poly_num: usize = 0;
            let mut max_coset_poly_num: usize = 0;
            for p in plans {
                log::debug!(
                    "## {}->{}",
                    coset_poly_num,
                    coset_poly_num + p.alloc_before.len()
                );
                coset_poly_num += p.alloc_before.len();
                if coset_poly_num > max_coset_poly_num {
                    max_coset_poly_num = coset_poly_num;
                }
                coset_polys.extend(p.alloc_before.iter());
                debug_assert_eq!(coset_polys.len(), coset_poly_num);
                log::debug!(
                    "## {}->{}",
                    coset_poly_num,
                    coset_poly_num - p.release_after.len()
                );
                coset_poly_num -= p.release_after.len();
                for r in &p.release_after {
                    let is_exist = coset_polys.remove(r);
                    assert!(is_exist);
                }
                debug_assert_eq!(coset_polys.len(), coset_poly_num);
            }
            debug_assert!(coset_polys.is_empty());
            (
                plans.iter().map(|p| p.alloc_before.len()).sum::<usize>(),
                max_coset_poly_num,
            )
        };
        log::debug!("simulate gates eval");
        let gates_info = simulate(&self.gates);
        log::debug!("simulate lookups eval");
        let lookup_info = simulate(&self.lookups);
        (
            gates_info.0 + lookup_info.0,
            std::cmp::max(gates_info.1, lookup_info.1),
        )
    }
}

impl<C: CurveAffine> Evaluator<C> {
    /// Creates a new evaluation structure
    pub fn new(cs: &ConstraintSystem<C::ScalarExt>) -> Self {
        let mem_plan = Self::analyse_coset_memory_plan(cs);
        let mut ev = Evaluator::default();

        // Custom gates
        for ev_group_id in 0..cs.evaluation_group_num() {
            let mut gev = GraphEvaluator::default();
            let mut parts = Vec::new();
            for gate in cs.gates.iter().filter(|g| g.ev_group_idx == ev_group_id) {
                parts.extend(
                    gate.polynomials()
                        .iter()
                        .map(|poly| gev.add_expression(poly)),
                );
            }
            gev.add_calculation(Calculation::Horner(
                ValueSource::PreviousValue(),
                parts,
                ValueSource::Y(),
            ));
            ev.custom_gates.push(gev);
        }

        // Lookups
        for ev_group_id in 0..cs.evaluation_group_num() {
            let mut evs = Vec::new();
            for lookup in cs.lookups.iter().filter(|l| l.ev_group_idx == ev_group_id) {
                let mut graph = GraphEvaluator::default();
                let mut evaluate_lc = |expressions: &Vec<Expression<_>>| {
                    let parts = expressions
                        .iter()
                        .map(|expr| graph.add_expression(expr))
                        .collect();
                    graph.add_calculation(Calculation::Horner(
                        ValueSource::Constant(0),
                        parts,
                        ValueSource::Theta(),
                    ))
                };

                // Input coset
                let compressed_input_coset = evaluate_lc(&lookup.input_expressions);
                // table coset
                let compressed_table_coset = evaluate_lc(&lookup.table_expressions);
                // z(\omega X) (a'(X) + \beta) (s'(X) + \gamma)
                let right_gamma = graph.add_calculation(Calculation::Add(
                    compressed_table_coset,
                    ValueSource::Gamma(),
                ));
                let lc = graph.add_calculation(Calculation::Add(
                    compressed_input_coset,
                    ValueSource::Beta(),
                ));
                graph.add_calculation(Calculation::Mul(lc, right_gamma));

                evs.push((graph, lookup.ev_group_idx));
            }
            ev.lookups.push(evs);
        }
        ev.mem_plan = mem_plan;
        ev
    }

    /// Evaluate h poly
    pub(in crate::plonk) fn evaluate_h(
        &self,
        pk: &ProvingKey<C>,
        advice_polys: &[&[Polynomial<C::ScalarExt, Coeff>]],
        instance_polys: &[&[Polynomial<C::ScalarExt, Coeff>]],
        y: C::ScalarExt,
        beta: C::ScalarExt,
        gamma: C::ScalarExt,
        theta: C::ScalarExt,
        lookups: &[Vec<lookup::prover::Committed<C>>],
        permutations: &[permutation::prover::Committed<C>],
    ) -> Polynomial<C::ScalarExt, ExtendedLagrangeCoeff> {
        let domain = &pk.vk.domain;
        log::debug!(
            "evaluate_h domain k: {} extended_k: {}",
            domain.k(),
            domain.extended_k()
        );
        let size = domain.extended_len();
        let rot_scale = 1 << (domain.extended_k() - domain.k());
        let fixed = &pk.fixed_cosets[..];
        let extended_omega = domain.get_extended_omega();
        let isize = size as i32;
        let one = C::ScalarExt::one();
        let l0 = &pk.l0;
        let l_last = &pk.l_last;
        let l_active_row = &pk.l_active_row;
        let p = &pk.vk.cs.permutation;

        // Calculate the advice and instance cosets
        let mut advice: Vec<Vec<Option<Polynomial<C::Scalar, ExtendedLagrangeCoeff>>>> =
            advice_polys
                .iter()
                .map(|advice_polys| advice_polys.iter().map(|_poly| None).collect())
                .collect();
        let mut instance: Vec<Vec<Option<Polynomial<C::Scalar, ExtendedLagrangeCoeff>>>> =
            instance_polys
                .iter()
                .map(|instance_polys| instance_polys.iter().map(|_poly| None).collect())
                .collect();

        let mut values = domain.empty_extended();

        // Core expression evaluations
        let num_threads = multicore::current_num_threads();
        for (cidx, (((advice, instance), lookups), permutation)) in advice
            .iter_mut()
            .zip(instance.iter_mut())
            .zip(lookups.iter())
            .zip(permutations.iter())
            .enumerate()
        {
            let alloc_cosets =
                |columns: &[circuit::Column<Any>],
                 advice: &mut Vec<Option<Polynomial<_, _>>>,
                 instance: &mut Vec<Option<Polynomial<_, _>>>| {
                    for c in columns {
                        match c.column_type {
                            Any::Advice => {
                                *advice.index_mut(c.index) = Some(
                                    domain.coeff_to_extended(advice_polys[cidx][c.index].clone()),
                                );
                            }
                            Any::Instance => {
                                *instance.index_mut(c.index) = Some(
                                    domain.coeff_to_extended(instance_polys[cidx][c.index].clone()),
                                );
                            }
                            Any::Fixed => {}
                        }
                    }
                };
            let release_cosets =
                |columns: &[circuit::Column<Any>],
                 advice: &mut Vec<Option<Polynomial<_, _>>>,
                 instance: &mut Vec<Option<Polynomial<_, _>>>| {
                    for c in columns {
                        match c.column_type {
                            Any::Advice => {
                                *advice.index_mut(c.index) = None;
                            }
                            Any::Instance => {
                                *instance.index_mut(c.index) = None;
                            }
                            Any::Fixed => {}
                        }
                    }
                };

            // Custom gates
            for (idx, gev) in self.custom_gates.iter().enumerate() {
                // first release unused coset ploys and calculated needed coset polys
                let plan = &self.mem_plan.gates[idx];
                alloc_cosets(&plan.alloc_before, advice, instance);
                let instance_readonly: &Vec<_> = instance;
                let advice_readonly: &Vec<_> = advice;

                multicore::scope(|scope| {
                    let chunk_size = (size + num_threads - 1) / num_threads;
                    for (thread_idx, values) in values.chunks_mut(chunk_size).enumerate() {
                        let start = thread_idx * chunk_size;
                        scope.spawn(move |_| {
                            let mut eval_data = gev.instance();
                            for (i, value) in values.iter_mut().enumerate() {
                                let idx = start + i;
                                *value = gev.evaluate(
                                    &mut eval_data,
                                    fixed,
                                    advice_readonly,
                                    instance_readonly,
                                    &beta,
                                    &gamma,
                                    &theta,
                                    &y,
                                    value,
                                    idx,
                                    rot_scale,
                                    isize,
                                );
                            }
                        });
                    }
                });
                release_cosets(&plan.release_after, advice, instance);
            }

            // Permutations
            let sets = &permutation.sets;
            if !sets.is_empty() {
                // TODO: optimize this
                alloc_cosets(&p.columns, advice, instance);

                let blinding_factors = pk.vk.cs.blinding_factors();
                let last_rotation = Rotation(-((blinding_factors + 1) as i32));
                let chunk_len = pk.vk.cs.degree() - 2;
                let delta_start = beta * &C::Scalar::ZETA;

                let first_set = sets.first().unwrap();
                let last_set = sets.last().unwrap();

                // Permutation constraints
                parallelize(&mut values, |values, start| {
                    let mut beta_term = extended_omega.pow_vartime(&[start as u64, 0, 0, 0]);
                    for (i, value) in values.iter_mut().enumerate() {
                        let idx = start + i;
                        let r_next = get_rotation_idx(idx, 1, rot_scale, isize);
                        let r_last = get_rotation_idx(idx, last_rotation.0, rot_scale, isize);

                        // Enforce only for the first set.
                        // l_0(X) * (1 - z_0(X)) = 0
                        *value = *value * y
                            + ((one - first_set.permutation_product_coset[idx]) * l0[idx]);
                        // Enforce only for the last set.
                        // l_last(X) * (z_l(X)^2 - z_l(X)) = 0
                        *value = *value * y
                            + ((last_set.permutation_product_coset[idx]
                                * last_set.permutation_product_coset[idx]
                                - last_set.permutation_product_coset[idx])
                                * l_last[idx]);
                        // Except for the first set, enforce.
                        // l_0(X) * (z_i(X) - z_{i-1}(\omega^(last) X)) = 0
                        for (set_idx, set) in sets.iter().enumerate() {
                            if set_idx != 0 {
                                *value = *value * y
                                    + ((set.permutation_product_coset[idx]
                                        - permutation.sets[set_idx - 1].permutation_product_coset
                                            [r_last])
                                        * l0[idx]);
                            }
                        }
                        // And for all the sets we enforce:
                        // (1 - (l_last(X) + l_blind(X))) * (
                        //   z_i(\omega X) \prod_j (p(X) + \beta s_j(X) + \gamma)
                        // - z_i(X) \prod_j (p(X) + \delta^j \beta X + \gamma)
                        // )
                        let mut current_delta = delta_start * beta_term;
                        for ((set, columns), cosets) in sets
                            .iter()
                            .zip(p.columns.chunks(chunk_len))
                            .zip(pk.permutation.cosets.chunks(chunk_len))
                        {
                            let mut left = set.permutation_product_coset[r_next];
                            for (values, permutation) in columns
                                .iter()
                                .map(|&column| match column.column_type() {
                                    Any::Fixed => &fixed[column.index()],
                                    Any::Advice => advice[column.index()].as_ref().unwrap(),
                                    Any::Instance => instance[column.index()].as_ref().unwrap(),
                                })
                                .zip(cosets.iter())
                            {
                                left *= values[idx] + beta * permutation[idx] + gamma;
                            }

                            let mut right = set.permutation_product_coset[idx];
                            for values in columns.iter().map(|&column| match column.column_type() {
                                Any::Fixed => &fixed[column.index()],
                                Any::Advice => advice[column.index()].as_ref().unwrap(),
                                Any::Instance => instance[column.index()].as_ref().unwrap(),
                            }) {
                                right *= values[idx] + current_delta + gamma;
                                current_delta *= &C::Scalar::DELTA;
                            }

                            *value = *value * y + ((left - right) * l_active_row[idx]);
                        }
                        beta_term *= &extended_omega;
                    }
                });

                release_cosets(&p.columns, advice, instance);
            }

            // Lookups
            for (ev_group, lookups) in &lookups
                .iter()
                .zip(self.lookups.iter().flatten())
                .group_by(|(_lookup, (_lookup_evaluator, ev_group))| ev_group)
            {
                // first release unused coset ploys and calculated needed coset polys
                let plan = &self.mem_plan.lookups[*ev_group];
                alloc_cosets(&plan.alloc_before, advice, instance);
                for (lookup, (lookup_evaluator, _ev_group)) in lookups {
                    // Polynomials required for this lookup.
                    // Calculated here so these only have to be kept in memory for the short time
                    // they are actually needed.
                    let product_coset = pk.vk.domain.coeff_to_extended(lookup.product_poly.clone());
                    let permuted_input_coset = pk
                        .vk
                        .domain
                        .coeff_to_extended(lookup.permuted_input_poly.clone());
                    let permuted_table_coset = pk
                        .vk
                        .domain
                        .coeff_to_extended(lookup.permuted_table_poly.clone());

                    // Lookup constraints
                    parallelize(&mut values, |values, start| {
                        let mut eval_data = lookup_evaluator.instance();
                        for (i, value) in values.iter_mut().enumerate() {
                            let idx = start + i;

                            let table_value = lookup_evaluator.evaluate(
                                &mut eval_data,
                                fixed,
                                advice,
                                instance,
                                &beta,
                                &gamma,
                                &theta,
                                &y,
                                &C::ScalarExt::zero(),
                                idx,
                                rot_scale,
                                isize,
                            );

                            let r_next = get_rotation_idx(idx, 1, rot_scale, isize);
                            let r_prev = get_rotation_idx(idx, -1, rot_scale, isize);

                            let a_minus_s = permuted_input_coset[idx] - permuted_table_coset[idx];
                            // l_0(X) * (1 - z(X)) = 0
                            *value = *value * y + ((one - product_coset[idx]) * l0[idx]);
                            // l_last(X) * (z(X)^2 - z(X)) = 0
                            *value = *value * y
                                + ((product_coset[idx] * product_coset[idx] - product_coset[idx])
                                    * l_last[idx]);
                            // (1 - (l_last(X) + l_blind(X))) * (
                            //   z(\omega X) (a'(X) + \beta) (s'(X) + \gamma)
                            //   - z(X) (\theta^{m-1} a_0(X) + ... + a_{m-1}(X) + \beta)
                            //          (\theta^{m-1} s_0(X) + ... + s_{m-1}(X) + \gamma)
                            // ) = 0
                            *value = *value * y
                                + ((product_coset[r_next]
                                    * (permuted_input_coset[idx] + beta)
                                    * (permuted_table_coset[idx] + gamma)
                                    - product_coset[idx] * table_value)
                                    * l_active_row[idx]);
                            // Check that the first values in the permuted input expression and permuted
                            // fixed expression are the same.
                            // l_0(X) * (a'(X) - s'(X)) = 0
                            *value = *value * y + (a_minus_s * l0[idx]);
                            // Check that each value in the permuted lookup input expression is either
                            // equal to the value above it, or the value at the same index in the
                            // permuted table expression.
                            // (1 - (l_last + l_blind)) * (a′(X) − s′(X))⋅(a′(X) − a′(\omega^{-1} X)) = 0
                            *value = *value * y
                                + (a_minus_s
                                    * (permuted_input_coset[idx] - permuted_input_coset[r_prev])
                                    * l_active_row[idx]);
                        }
                    });
                }
                release_cosets(&plan.release_after, advice, instance);
            }
        }
        values
    }

    fn analyse_coset_memory_plan(cs: &ConstraintSystem<C::ScalarExt>) -> MemPlan {
        use itertools::Itertools;
        let mut mem_plan = MemPlan::default();

        // step1: analyse gates
        let mut col_and_ev_group_of_gates = Vec::new();
        for gate in cs.gates.iter() {
            for poly in gate.polynomials() {
                let cols_used = Self::collect_columns_needed_by_expression(poly);
                for c in cols_used {
                    col_and_ev_group_of_gates.push((c, gate.ev_group_idx));
                }
            }
        }
        mem_plan.gates = Self::col_and_ev_group_to_mem_plan(cs, &col_and_ev_group_of_gates);

        // step2: analyse lookups
        let mut col_and_ev_group_of_lookups: Vec<(circuit::Column<Any>, usize)> = Vec::new();
        for lookup in &cs.lookups {
            let mut s: BTreeSet<circuit::Column<Any>> = Default::default();
            for e in lookup
                .input_expressions
                .iter()
                .chain(lookup.table_expressions.iter())
            {
                let cols_used = Self::collect_columns_needed_by_expression(e);
                s.extend(cols_used.iter())
            }
            for c in &s {
                col_and_ev_group_of_lookups.push((*c, lookup.ev_group_idx));
            }
        }
        mem_plan.lookups = Self::col_and_ev_group_to_mem_plan(cs, &col_and_ev_group_of_lookups);

        // print debug info of this plan
        {
            let (fft_num, max_coset_poly_num) = mem_plan.simulate_execution();
            let permutation_column_num = cs
                .permutation
                .columns
                .iter()
                .filter(|c| matches!(c.column_type, Any::Advice | Any::Instance))
                .count();

            log::debug!(
                "with coset mem opt: max coset poly num: {}",
                max_coset_poly_num
            );
            log::debug!(
                "with coset mem opt: total fft: {}",
                fft_num + permutation_column_num
            );

            log::debug!(
                "without coset mem opt: max coset poly num: {}",
                cs.num_advice_columns + cs.num_instance_columns
            );
            log::debug!(
                "without coset mem opt: total fft: {}",
                cs.num_advice_columns + cs.num_instance_columns
            );
        }

        mem_plan
    }

    fn collect_columns_needed_by_expression(
        e: &Expression<C::ScalarExt>,
    ) -> Vec<circuit::Column<Any>> {
        e.evaluate(
            &|_constant| Vec::new(),
            &|_selector| Vec::new(),
            &|_fixed| Vec::new(),
            &|advice| {
                vec![circuit::Column {
                    index: advice.column_index(),
                    column_type: Any::Advice,
                }]
            },
            &|instance| {
                vec![circuit::Column {
                    index: instance.column_index(),
                    column_type: Any::Instance,
                }]
            },
            &|negate| negate,
            &|sum1, sum2| sum1.into_iter().chain(sum2.into_iter()).unique().collect(),
            &|prod1, prod2| {
                prod1
                    .into_iter()
                    .chain(prod2.into_iter())
                    .unique()
                    .collect()
            },
            &|exp, _scalar| exp,
        )
    }

    fn col_and_ev_group_to_mem_plan(
        cs: &ConstraintSystem<C::ScalarExt>,
        col_and_group: &[(circuit::Column<Any>, usize)],
    ) -> Vec<MemoryPlanOfEvalutionGroup> {
        #[derive(Debug)]
        struct Lifetime {
            first_group: usize,
            last_group: usize,
        }
        let mut col_lifetime: Vec<(circuit::Column<Any>, Lifetime)> = Default::default();
        for (c, g) in col_and_group {
            let mut found: bool = false;
            let c = *c;
            for item in col_lifetime.iter_mut() {
                if item.0 == c {
                    item.1.last_group = *g;
                    found = true;
                }
            }
            if !found {
                col_lifetime.push((
                    c,
                    Lifetime {
                        first_group: *g,
                        last_group: *g,
                    },
                ))
            }
        }

        let mut plans: Vec<MemoryPlanOfEvalutionGroup> = Vec::new();

        for t in 0..cs.evaluation_group_num() {
            let mut p = MemoryPlanOfEvalutionGroup::default();
            for (c, l) in &col_lifetime {
                if l.first_group == t {
                    p.alloc_before.push(*c);
                }
                if l.last_group == t {
                    p.release_after.push(*c);
                }
            }
            plans.push(p);
        }

        plans
    }
}

impl<C: CurveAffine> Default for GraphEvaluator<C> {
    fn default() -> Self {
        Self {
            // Fixed positions to allow easy access
            constants: vec![
                C::ScalarExt::zero(),
                C::ScalarExt::one(),
                C::ScalarExt::from(2u64),
            ],
            rotations: Vec::new(),
            calculations: Vec::new(),
            num_intermediates: 0,
        }
    }
}

impl<C: CurveAffine> GraphEvaluator<C> {
    /// Adds a rotation
    fn add_rotation(&mut self, rotation: &Rotation) -> usize {
        let position = self.rotations.iter().position(|&c| c == rotation.0);
        match position {
            Some(pos) => pos,
            None => {
                self.rotations.push(rotation.0);
                self.rotations.len() - 1
            }
        }
    }

    /// Adds a constant
    fn add_constant(&mut self, constant: &C::ScalarExt) -> ValueSource {
        let position = self.constants.iter().position(|&c| c == *constant);
        ValueSource::Constant(match position {
            Some(pos) => pos,
            None => {
                self.constants.push(*constant);
                self.constants.len() - 1
            }
        })
    }

    /// Adds a calculation.
    /// Currently does the simplest thing possible: just stores the
    /// resulting value so the result can be reused  when that calculation
    /// is done multiple times.
    fn add_calculation(&mut self, calculation: Calculation) -> ValueSource {
        let existing_calculation = self
            .calculations
            .iter()
            .find(|c| c.calculation == calculation);
        match existing_calculation {
            Some(existing_calculation) => ValueSource::Intermediate(existing_calculation.target),
            None => {
                let target = self.num_intermediates;
                self.calculations.push(CalculationInfo {
                    calculation,
                    target,
                });
                self.num_intermediates += 1;
                ValueSource::Intermediate(target)
            }
        }
    }

    /// Generates an optimized evaluation for the expression
    fn add_expression(&mut self, expr: &Expression<C::ScalarExt>) -> ValueSource {
        match expr {
            Expression::Constant(scalar) => self.add_constant(scalar),
            Expression::Selector(_selector) => unreachable!(),
            Expression::Fixed(query) => {
                let rot_idx = self.add_rotation(&query.rotation);
                self.add_calculation(Calculation::Store(ValueSource::Fixed(
                    query.column_index,
                    rot_idx,
                )))
            }
            Expression::Advice(query) => {
                let rot_idx = self.add_rotation(&query.rotation);
                self.add_calculation(Calculation::Store(ValueSource::Advice(
                    query.column_index,
                    rot_idx,
                )))
            }
            Expression::Instance(query) => {
                let rot_idx = self.add_rotation(&query.rotation);
                self.add_calculation(Calculation::Store(ValueSource::Instance(
                    query.column_index,
                    rot_idx,
                )))
            }
            Expression::Negated(a) => match **a {
                Expression::Constant(scalar) => self.add_constant(&-scalar),
                _ => {
                    let result_a = self.add_expression(a);
                    match result_a {
                        ValueSource::Constant(0) => result_a,
                        _ => self.add_calculation(Calculation::Negate(result_a)),
                    }
                }
            },
            Expression::Sum(a, b) => {
                // Undo subtraction stored as a + (-b) in expressions
                match &**b {
                    Expression::Negated(b_int) => {
                        let result_a = self.add_expression(a);
                        let result_b = self.add_expression(b_int);
                        if result_a == ValueSource::Constant(0) {
                            self.add_calculation(Calculation::Negate(result_b))
                        } else if result_b == ValueSource::Constant(0) {
                            result_a
                        } else {
                            self.add_calculation(Calculation::Sub(result_a, result_b))
                        }
                    }
                    _ => {
                        let result_a = self.add_expression(a);
                        let result_b = self.add_expression(b);
                        if result_a == ValueSource::Constant(0) {
                            result_b
                        } else if result_b == ValueSource::Constant(0) {
                            result_a
                        } else if result_a <= result_b {
                            self.add_calculation(Calculation::Add(result_a, result_b))
                        } else {
                            self.add_calculation(Calculation::Add(result_b, result_a))
                        }
                    }
                }
            }
            Expression::Product(a, b) => {
                let result_a = self.add_expression(a);
                let result_b = self.add_expression(b);
                if result_a == ValueSource::Constant(0) || result_b == ValueSource::Constant(0) {
                    ValueSource::Constant(0)
                } else if result_a == ValueSource::Constant(1) {
                    result_b
                } else if result_b == ValueSource::Constant(1) {
                    result_a
                } else if result_a == ValueSource::Constant(2) {
                    self.add_calculation(Calculation::Double(result_b))
                } else if result_b == ValueSource::Constant(2) {
                    self.add_calculation(Calculation::Double(result_a))
                } else if result_a == result_b {
                    self.add_calculation(Calculation::Square(result_a))
                } else if result_a <= result_b {
                    self.add_calculation(Calculation::Mul(result_a, result_b))
                } else {
                    self.add_calculation(Calculation::Mul(result_b, result_a))
                }
            }
            Expression::Scaled(a, f) => {
                if *f == C::ScalarExt::zero() {
                    ValueSource::Constant(0)
                } else if *f == C::ScalarExt::one() {
                    self.add_expression(a)
                } else {
                    let cst = self.add_constant(f);
                    let result_a = self.add_expression(a);
                    self.add_calculation(Calculation::Mul(result_a, cst))
                }
            }
        }
    }

    /// Creates a new evaluation structure
    pub fn instance(&self) -> EvaluationData<C> {
        EvaluationData {
            intermediates: vec![C::ScalarExt::zero(); self.num_intermediates],
            rotations: vec![0usize; self.rotations.len()],
        }
    }

    pub fn evaluate<B: Basis>(
        &self,
        data: &mut EvaluationData<C>,
        fixed: &[Polynomial<C::ScalarExt, B>],
        advice: &[Option<Polynomial<C::ScalarExt, B>>],
        instance: &[Option<Polynomial<C::ScalarExt, B>>],
        beta: &C::ScalarExt,
        gamma: &C::ScalarExt,
        theta: &C::ScalarExt,
        y: &C::ScalarExt,
        previous_value: &C::ScalarExt,
        idx: usize,
        rot_scale: i32,
        isize: i32,
    ) -> C::ScalarExt {
        // All rotation index values
        for (rot_idx, rot) in self.rotations.iter().enumerate() {
            data.rotations[rot_idx] = get_rotation_idx(idx, *rot, rot_scale, isize);
        }

        // All calculations, with cached intermediate results
        for calc in self.calculations.iter() {
            data.intermediates[calc.target] = calc.calculation.evaluate(
                &data.rotations,
                &self.constants,
                &data.intermediates,
                fixed,
                advice,
                instance,
                beta,
                gamma,
                theta,
                y,
                previous_value,
            );
        }

        // Return the result of the last calculation (if any)
        if let Some(calc) = self.calculations.last() {
            data.intermediates[calc.target]
        } else {
            C::ScalarExt::zero()
        }
    }
}

/// Simple evaluation of an expression
pub fn evaluate<F: FieldExt, B: Basis>(
    expression: &Expression<F>,
    size: usize,
    rot_scale: i32,
    fixed: &[Polynomial<F, B>],
    advice: &[Polynomial<F, B>],
    instance: &[Polynomial<F, B>],
) -> Vec<F> {
    let mut values = vec![F::zero(); size];
    let isize = size as i32;
    parallelize(&mut values, |values, start| {
        for (i, value) in values.iter_mut().enumerate() {
            let idx = start + i;
            *value = expression.evaluate(
                &|scalar| scalar,
                &|_| panic!("virtual selectors are removed during optimization"),
                &|query| {
                    fixed[query.column_index]
                        [get_rotation_idx(idx, query.rotation.0, rot_scale, isize)]
                },
                &|query| {
                    advice[query.column_index]
                        [get_rotation_idx(idx, query.rotation.0, rot_scale, isize)]
                },
                &|query| {
                    instance[query.column_index]
                        [get_rotation_idx(idx, query.rotation.0, rot_scale, isize)]
                },
                &|a| -a,
                &|a, b| a + &b,
                &|a, b| a * b,
                &|a, scalar| a * scalar,
            );
        }
    });
    values
}
