use ff::Field;
use group::Curve;
use rand_core::RngCore;
use std::iter;

use super::{
    vanishing, ChallengeBeta, ChallengeGamma, ChallengeTheta, ChallengeX, ChallengeY, Error,
    VerifyingKey,
};
use crate::arithmetic::{compute_inner_product, CurveAffine, FieldExt};
use crate::poly::commitment::{CommitmentScheme, Verifier};
use crate::poly::VerificationStrategy;
use crate::poly::{
    commitment::{Blind, Params, MSM},
    Guard, VerifierQuery,
};
use crate::transcript::{self, read_n_points, read_n_scalars, EncodedChallenge, TranscriptRead};

#[cfg(feature = "batch")]
mod batch;
#[cfg(feature = "batch")]
pub use batch::BatchVerifier;

use crate::poly::commitment::ParamsVerifier;

/// Returns a boolean indicating whether or not the proof is valid
pub fn verify_proof<
    'params,
    Scheme: CommitmentScheme,
    V: Verifier<'params, Scheme>,
    E: EncodedChallenge<Scheme::Curve>,
    T: TranscriptRead<Scheme::Curve, E>,
    Strategy: VerificationStrategy<'params, Scheme, V>,
>(
    params: &'params Scheme::ParamsVerifier,
    vk: &VerifyingKey<Scheme::Curve>,
    strategy: Strategy,
    instances: &[&[&[Scheme::Scalar]]],
    transcript: &mut T,
) -> Result<Strategy::Output, Error> {
    // Check that instances matches the expected number of instance columns
    for instances in instances.iter() {
        if instances.len() != vk.cs.num_instance_columns {
            return Err(Error::InvalidInstances);
        }
    }

    let instance_commitments = if V::QUERY_INSTANCE {
        instances
            .iter()
            .map(|instance| {
                instance
                    .iter()
                    .map(|instance| {
                        if instance.len() > params.n() as usize - (vk.cs.blinding_factors() + 1) {
                            return Err(Error::InstanceTooLarge);
                        }
                        let mut poly = instance.to_vec();
                        poly.resize(params.n() as usize, Scheme::Scalar::zero());
                        let poly = vk.domain.lagrange_from_vec(poly);

                        Ok(params.commit_lagrange(&poly, Blind::default()).to_affine())
                    })
                    .collect::<Result<Vec<_>, _>>()
            })
            .collect::<Result<Vec<_>, _>>()?
    } else {
        vec![vec![]; instances.len()]
    };

    println!("instance_commitments");
    for instance_commitments in instance_commitments.iter() {
        for instance_commitment in instance_commitments.iter() {
            println!("{:?}", instance_commitment);
        }
    }

    let num_proofs = instance_commitments.len();

    // Hash verification key into transcript
    vk.hash_into(transcript)?;

    if V::QUERY_INSTANCE {
        for instance_commitments in instance_commitments.iter() {
            // Hash the instance (external) commitments into the transcript
            for commitment in instance_commitments {
                transcript.common_point(*commitment)?
            }
        }
    } else {
        for instance in instances.iter() {
            for instance in instance.iter() {
                for value in instance.iter() {
                    transcript.common_scalar(*value)?;
                }
            }
        }
    }

    // Hash the prover's advice commitments into the transcript and squeeze challenges
    let (advice_commitments, challenges) = {
        let mut advice_commitments =
            vec![vec![Scheme::Curve::default(); vk.cs.num_advice_columns]; num_proofs];
        let mut challenges = vec![Scheme::Scalar::zero(); vk.cs.num_challenges];

        for current_phase in vk.cs.phases() {
            for advice_commitments in advice_commitments.iter_mut() {
                for (phase, commitment) in vk
                    .cs
                    .advice_column_phase
                    .iter()
                    .zip(advice_commitments.iter_mut())
                {
                    if current_phase == *phase {
                        *commitment = transcript.read_point()?;
                    }
                }
            }
            for (phase, challenge) in vk.cs.challenge_phase.iter().zip(challenges.iter_mut()) {
                if current_phase == *phase {
                    *challenge = *transcript.squeeze_challenge_scalar::<()>();
                }
            }
        }

        (advice_commitments, challenges)
    };

    println!("expected_advice_commitments_vec");
    for advice_commitment in advice_commitments.iter() {
        for advice_commitment in advice_commitment.iter() {
            // println!("{:?}", advice_commitment);
        }
    }
    println!("challenges");
    for challenge in challenges.iter() {
        println!("{:?}", challenge);
    }

    // Sample theta challenge for keeping lookup columns linearly independent
    let theta: ChallengeTheta<_> = transcript.squeeze_challenge_scalar();
    println!("theta: {:?}", *theta);

    let lookups_permuted = (0..num_proofs)
        .map(|_| -> Result<Vec<_>, _> {
            // Hash each lookup permuted commitment
            vk.cs
                .lookups
                .iter()
                .map(|argument| argument.read_permuted_commitments(transcript))
                .collect::<Result<Vec<_>, _>>()
        })
        .collect::<Result<Vec<_>, _>>()?;

    println!("lookup_permuted_commitments_vec");
    for lookup_permuted in lookups_permuted.iter() {
        println!("lookup_permuted_commitment");
        for lookup in lookup_permuted.iter() {
            println!("input: {:?}", lookup.permuted_input_commitment);
            println!("table: {:?}", lookup.permuted_table_commitment);
        }
    }

    // Sample beta challenge
    let beta: ChallengeBeta<_> = transcript.squeeze_challenge_scalar();

    // Sample gamma challenge
    let gamma: ChallengeGamma<_> = transcript.squeeze_challenge_scalar();

    println!("beta: {:?}", *beta);
    println!("gamma: {:?}", *gamma);

    let permutations_committed = (0..num_proofs)
        .map(|_| {
            // Hash each permutation product commitment
            vk.cs.permutation.read_product_commitments(vk, transcript)
        })
        .collect::<Result<Vec<_>, _>>()?;

    println!("expected_permutation_product_commitments_vec");
    for permutations_committed in permutations_committed.iter() {
        println!("permutation_committed");
        for permutation_product_commitment in permutations_committed
            .permutation_product_commitments
            .iter()
        {
            println!("product_commitment: {:?}", *permutation_product_commitment);
        }
    }

    let lookups_committed = lookups_permuted
        .into_iter()
        .map(|lookups| {
            // Hash each lookup product commitment
            lookups
                .into_iter()
                .map(|lookup| lookup.read_product_commitment(transcript))
                .collect::<Result<Vec<_>, _>>()
        })
        .collect::<Result<Vec<_>, _>>()?;

    println!("lookup_product_commitments_vec");
    for lookup_committed in lookups_committed.iter() {
        println!("lookup_committed");
        for lookup_committed in lookup_committed.iter() {
            println!(
                "product_commitment: {:?}",
                lookup_committed.product_commitment
            );
        }
    }

    let vanishing = vanishing::Argument::read_commitments_before_y(transcript)?;

    println!("expected_vanishing_random_poly_commitment: {:?}", vanishing.random_poly_commitment);

    // Sample y challenge, which keeps the gates linearly independent.
    let y: ChallengeY<_> = transcript.squeeze_challenge_scalar();
    println!("expected_y: {:?}", *y);

    let vanishing = vanishing.read_commitments_after_y(vk, transcript)?;

    println!("expected_vanishing_h_poly_commitments");
    for h_commitment in vanishing.h_commitments.iter() {
        println!("h_commitment: {:?}", h_commitment);
    }

    // Sample x challenge, which is used to ensure the circuit is
    // satisfied with high probability.
    let x: ChallengeX<_> = transcript.squeeze_challenge_scalar();
    println!("expected_x: {:?}", *x);
    let instance_evals = if V::QUERY_INSTANCE {
        (0..num_proofs)
            .map(|_| -> Result<Vec<_>, _> {
                read_n_scalars(transcript, vk.cs.instance_queries.len())
            })
            .collect::<Result<Vec<_>, _>>()?
    } else {
        let xn = x.pow(&[params.n() as u64, 0, 0, 0]);
        let (min_rotation, max_rotation) =
            vk.cs
                .instance_queries
                .iter()
                .fold((0, 0), |(min, max), (_, rotation)| {
                    if rotation.0 < min {
                        (rotation.0, max)
                    } else if rotation.0 > max {
                        (min, rotation.0)
                    } else {
                        (min, max)
                    }
                });
        let max_instance_len = instances
            .iter()
            .flat_map(|instance| instance.iter().map(|instance| instance.len()))
            .max_by(Ord::cmp)
            .unwrap_or_default();
        let l_i_s = &vk.domain.l_i_range(
            *x,
            xn,
            -max_rotation..max_instance_len as i32 + min_rotation.abs(),
        );
        instances
            .iter()
            .map(|instances| {
                vk.cs
                    .instance_queries
                    .iter()
                    .map(|(column, rotation)| {
                        let instances = instances[column.index()];
                        let offset = (max_rotation - rotation.0) as usize;
                        compute_inner_product(instances, &l_i_s[offset..offset + instances.len()])
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>()
    };

    let advice_evals = (0..num_proofs)
        .map(|_| -> Result<Vec<_>, _> { read_n_scalars(transcript, vk.cs.advice_queries.len()) })
        .collect::<Result<Vec<_>, _>>()?;

    println!("expected_advice_evals_vec");
    for advice_evals in advice_evals.iter() {
        println!("advice_eval");
        for advice_eval in advice_evals.iter() {
            println!("{:?}", advice_eval);
        }
    }

    let fixed_evals = read_n_scalars(transcript, vk.cs.fixed_queries.len())?;

    println!("expected_fixed_evals_vec");
    for fixed_eval in fixed_evals.iter() {
        println!("{:?}", fixed_eval);
    }

    let vanishing = vanishing.evaluate_after_x(transcript)?;

    println!("expected_random_vanishing_eval: {:?}", vanishing.random_eval);

    let permutations_common = vk.permutation.evaluate(transcript)?;

    println!("expected_common_permutation_evals");
    for permutation_eval in permutations_common.permutation_evals.iter() {
        println!("permutation_eval: {:?}", permutation_eval);
    }

    let permutations_evaluated = permutations_committed
        .into_iter()
        .map(|permutation| permutation.evaluate(transcript))
        .collect::<Result<Vec<_>, _>>()?;

    println!("expected_permutation_evals_vec");
    for permutations_evaluated in permutations_evaluated.iter() {
        for permutation_eval_set in permutations_evaluated.sets.iter() {
            println!("eval: {:?}", permutation_eval_set.permutation_product_eval);
            println!(
                "next_eval: {:?}",
                permutation_eval_set.permutation_product_next_eval
            );
            if let Some(last) = permutation_eval_set.permutation_product_last_eval {
                println!("last: {:?}", last);
            } else {
                println!("empty");
            }
        }
    }

    println!("expected_lookup_evals_vec");
    let lookups_evaluated = lookups_committed
        .into_iter()
        .map(|lookups| -> Result<Vec<_>, _> {
            lookups
                .into_iter()
                .map(|lookup| lookup.evaluate(transcript))
                .collect::<Result<Vec<_>, _>>()
        })
        .collect::<Result<Vec<_>, _>>()?;

    println!("lookup_eval(eval, next, last)");
    for lookup_evaluated in lookups_evaluated.iter() {
        for evaluated in lookup_evaluated.iter() {
            println!("product_eval: {:?}", evaluated.product_eval);
            println!("product_next_eval: {:?}", evaluated.product_next_eval);
            println!("permuted_input_eval: {:?}", evaluated.permuted_input_eval);
            println!(
                "permuted_input_inv_eval: {:?}",
                evaluated.permuted_input_inv_eval
            );
            println!("permuted_table_eval: {:?}", evaluated.permuted_table_eval);
        }
    }

    // This check ensures the circuit is satisfied so long as the polynomial
    // commitments open to the correct values.
    let vanishing = {
        // x^n
        let xn = x.pow(&[params.n() as u64, 0, 0, 0]);

        let blinding_factors = vk.cs.blinding_factors();
        let l_evals = vk
            .domain
            .l_i_range(*x, xn, (-((blinding_factors + 1) as i32))..=0);
        assert_eq!(l_evals.len(), 2 + blinding_factors);
        let l_last = l_evals[0];
        let l_blind: Scheme::Scalar = l_evals[1..(1 + blinding_factors)]
            .iter()
            .fold(Scheme::Scalar::zero(), |acc, eval| acc + eval);
        let l_0 = l_evals[1 + blinding_factors];

        // Compute the expected value of h(x)
        let expressions = advice_evals
            .iter()
            .zip(instance_evals.iter())
            .zip(permutations_evaluated.iter())
            .zip(lookups_evaluated.iter())
            .flat_map(|(((advice_evals, instance_evals), permutation), lookups)| {
                let challenges = &challenges;
                let fixed_evals = &fixed_evals;
                std::iter::empty()
                    // Evaluate the circuit using the custom gates provided
                    .chain(vk.cs.gates.iter().flat_map(move |gate| {
                        gate.polynomials().iter().map(move |poly| {
                            poly.evaluate(
                                &|scalar| scalar,
                                &|_| panic!("virtual selectors are removed during optimization"),
                                &|query| fixed_evals[query.index],
                                &|query| advice_evals[query.index],
                                &|query| instance_evals[query.index],
                                &|challenge| challenges[challenge.index()],
                                &|a| -a,
                                &|a, b| a + &b,
                                &|a, b| a * &b,
                                &|a, scalar| a * &scalar,
                            )
                        })
                    }))
                    .chain(permutation.expressions(
                        vk,
                        &vk.cs.permutation,
                        &permutations_common,
                        advice_evals,
                        fixed_evals,
                        instance_evals,
                        l_0,
                        l_last,
                        l_blind,
                        beta,
                        gamma,
                        x,
                    ))
                    .chain(
                        lookups
                            .iter()
                            .zip(vk.cs.lookups.iter())
                            .flat_map(move |(p, argument)| {
                                p.expressions(
                                    l_0,
                                    l_last,
                                    l_blind,
                                    argument,
                                    theta,
                                    beta,
                                    gamma,
                                    advice_evals,
                                    fixed_evals,
                                    instance_evals,
                                    challenges,
                                )
                            })
                            .into_iter(),
                    )
            });

        vanishing.verify(params, expressions, y, xn)
    };

    let queries = instance_commitments
        .iter()
        .zip(instance_evals.iter())
        .zip(advice_commitments.iter())
        .zip(advice_evals.iter())
        .zip(permutations_evaluated.iter())
        .zip(lookups_evaluated.iter())
        .flat_map(
            |(
                (
                    (((instance_commitments, instance_evals), advice_commitments), advice_evals),
                    permutation,
                ),
                lookups,
            )| {
                iter::empty()
                    .chain(
                        V::QUERY_INSTANCE
                            .then_some(vk.cs.instance_queries.iter().enumerate().map(
                                move |(query_index, &(column, at))| {
                                    VerifierQuery::new_commitment(
                                        &instance_commitments[column.index()],
                                        vk.domain.rotate_omega(*x, at),
                                        instance_evals[query_index],
                                    )
                                },
                            ))
                            .into_iter()
                            .flatten(),
                    )
                    .chain(vk.cs.advice_queries.iter().enumerate().map(
                        move |(query_index, &(column, at))| {
                            VerifierQuery::new_commitment(
                                &advice_commitments[column.index()],
                                vk.domain.rotate_omega(*x, at),
                                advice_evals[query_index],
                            )
                        },
                    ))
                    .chain(permutation.queries(vk, x))
                    .chain(
                        lookups
                            .iter()
                            .flat_map(move |p| p.queries(vk, x))
                            .into_iter(),
                    )
            },
        )
        .chain(
            vk.cs
                .fixed_queries
                .iter()
                .enumerate()
                .map(|(query_index, &(column, at))| {
                    VerifierQuery::new_commitment(
                        &vk.fixed_commitments[column.index()],
                        vk.domain.rotate_omega(*x, at),
                        fixed_evals[query_index],
                    )
                }),
        )
        .chain(permutations_common.queries(&vk.permutation, x))
        .chain(vanishing.queries(x));

    // We are now convinced the circuit is satisfied so long as the
    // polynomial commitments open to the correct values.

    let verifier = V::new(params);
    strategy.process(|msm| {
        verifier
            .verify_proof(transcript, queries, msm)
            .map_err(|_| Error::Opening)
    })
}
