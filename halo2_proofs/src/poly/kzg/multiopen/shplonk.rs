mod prover;
mod verifier;

pub use prover::ProverSHPLONK;
pub use verifier::VerifierSHPLONK;

use crate::{
    arithmetic::{eval_polynomial, lagrange_interpolate, CurveAffine, FieldExt},
    poly::{query::Query, Coeff, Polynomial},
    transcript::ChallengeScalar,
};

use std::{
    collections::{btree_map::Entry, BTreeMap, BTreeSet},
    marker::PhantomData,
};

#[derive(Clone, Copy, Debug)]
struct U {}
type ChallengeU<F> = ChallengeScalar<F, U>;

#[derive(Clone, Copy, Debug)]
struct V {}
type ChallengeV<F> = ChallengeScalar<F, V>;

#[derive(Clone, Copy, Debug)]
struct Y {}
type ChallengeY<F> = ChallengeScalar<F, Y>;

#[derive(Debug, Clone, PartialEq)]
struct Commitment<F: FieldExt, T: PartialEq + Clone>((T, Vec<F>));

impl<F: FieldExt, T: PartialEq + Clone> Commitment<F, T> {
    fn get(&self) -> T {
        self.0 .0.clone()
    }

    fn evals(&self) -> Vec<F> {
        self.0 .1.clone()
    }
}

#[derive(Debug, Clone, PartialEq)]
struct RotationSet<F: FieldExt, T: PartialEq + Clone> {
    commitments: Vec<Commitment<F, T>>,
    points: Vec<F>,
}

#[derive(Debug, PartialEq)]
struct IntermediateSets<F: FieldExt, Q: Query<F>> {
    rotation_sets: Vec<RotationSet<F, Q::Commitment>>,
    super_point_set: Vec<F>,
}

fn construct_intermediate_sets<F: FieldExt, I, Q: Query<F, Eval = F>>(
    queries: I,
) -> IntermediateSets<F, Q>
where
    I: IntoIterator<Item = Q> + Clone,
{
    let queries = queries.into_iter().collect::<Vec<_>>();

    // Find evaluation of a commitment at a rotation
    let get_eval = |commitment: Q::Commitment, rotation: F| -> F {
        queries
            .iter()
            .find(|query| query.get_commitment() == commitment && query.get_point() == rotation)
            .unwrap()
            .get_eval()
    };

    // Order points according to their rotation
    let mut rotation_point_map = BTreeMap::new();
    for query in queries.clone() {
        let point = rotation_point_map
            .entry(query.get_point())
            .or_insert_with(|| query.get_point());

        // Assert rotation point matching consistency
        assert_eq!(*point, query.get_point());
    }
    // All points appear in queries
    let super_point_set: Vec<F> = rotation_point_map.values().cloned().collect();

    // Collect rotation sets for each commitment
    // Example elements in the vector:
    // (C_0, {r_5}),
    // (C_1, {r_1, r_2, r_3}),
    // (C_2, {r_2, r_3, r_4}),
    // (C_3, {r_2, r_3, r_4}),
    // ...
    let mut commitment_rotation_set_map: Vec<(Q::Commitment, BTreeSet<F>)> = vec![];
    for query in queries.clone() {
        let rotation = query.get_point();
        if let Some(pos) = commitment_rotation_set_map
            .iter()
            .position(|(commitment, _)| *commitment == query.get_commitment())
        {
            let (_, rotation_set) = &mut commitment_rotation_set_map[pos];
            rotation_set.insert(rotation);
        } else {
            let rotation_set = BTreeSet::from([rotation]);
            commitment_rotation_set_map.push((query.get_commitment(), rotation_set));
        };
    }

    // Flatten rotation sets and collect commitments that opens against each commitment set
    // Example elements in the vector:
    // {r_5}: [C_0],
    // {r_1, r_2, r_3} : [C_1]
    // {r_2, r_3, r_4} : [C_2, C_3],
    // ...
    let mut rotation_set_commitment_map = BTreeMap::<BTreeSet<_>, Vec<Q::Commitment>>::new();
    for (commitment, rotation_set) in commitment_rotation_set_map.iter() {
        let commitments = rotation_set_commitment_map
            .entry(rotation_set.clone())
            .or_insert_with(Vec::new);
        if !commitments.contains(commitment) {
            commitments.push(commitment.clone());
        }
    }

    let rotation_sets = rotation_set_commitment_map
        .into_iter()
        .map(|(rotation_set, commitments)| {
            let rotations: Vec<F> = rotation_set.iter().cloned().collect();

            let commitments: Vec<Commitment<F, Q::Commitment>> = commitments
                .iter()
                .map(|commitment| {
                    let evals: Vec<F> = rotations
                        .iter()
                        .map(|rotation| get_eval(commitment.clone(), *rotation))
                        .collect();
                    Commitment((commitment.clone(), evals))
                })
                .collect();

            RotationSet {
                commitments,
                points: rotations
                    .iter()
                    .map(|rotation| *rotation_point_map.get(rotation).unwrap())
                    .collect(),
            }
        })
        .collect::<Vec<RotationSet<_, _>>>();

    IntermediateSets {
        rotation_sets,
        super_point_set,
    }
}
