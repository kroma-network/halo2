use super::{
    construct_intermediate_sets, ChallengeU, ChallengeV, ChallengeY, Commitment, Query, RotationSet,
};
use crate::arithmetic::{
    eval_polynomial, evaluate_vanishing_polynomial, kate_division, lagrange_interpolate,
    CurveAffine, FieldExt,
};

use crate::poly::commitment::{Blind, ParamsProver, Prover};
use crate::poly::kzg::commitment::{KZGCommitmentScheme, ParamsKZG};
use crate::poly::query::{PolynomialPointer, ProverQuery};
use crate::poly::Rotation;
use crate::poly::{commitment::Params, Coeff, Polynomial};
use crate::transcript::{EncodedChallenge, TranscriptWrite};

use ff::Field;
use group::Curve;
use halo2curves::pairing::Engine;
use rand_core::RngCore;
use std::fmt::Debug;
use std::io::{self, Write};
use std::marker::PhantomData;
use std::ops::MulAssign;

fn div_by_vanishing<F: FieldExt>(poly: Polynomial<F, Coeff>, roots: &[F]) -> Vec<F> {
    let poly = roots
        .iter()
        .fold(poly.values, |poly, point| kate_division(&poly, *point));

    poly
}

struct CommitmentExtension<'a, C: CurveAffine> {
    commitment: Commitment<C::Scalar, PolynomialPointer<'a, C>>,
    low_degree_equivalent: Polynomial<C::Scalar, Coeff>,
}

impl<'a, C: CurveAffine> Commitment<C::Scalar, PolynomialPointer<'a, C>> {
    fn extend(&self, n: u64, points: Vec<C::Scalar>) -> CommitmentExtension<'a, C> {
        let mut poly = lagrange_interpolate(&points[..], &self.evals()[..]);
        poly.resize(n as usize, C::Scalar::zero());

        let low_degree_equivalent = Polynomial {
            values: poly,
            _marker: PhantomData,
        };

        CommitmentExtension {
            commitment: self.clone(),
            low_degree_equivalent,
        }
    }
}

impl<'a, C: CurveAffine> CommitmentExtension<'a, C> {
    fn linearisation_contribution(&self, u: C::Scalar) -> Polynomial<C::Scalar, Coeff> {
        let p_x = self.commitment.get().poly;
        let r_eval = eval_polynomial(&self.low_degree_equivalent.values[..], u);
        p_x - r_eval
    }

    fn quotient_contribution(&self) -> Polynomial<C::Scalar, Coeff> {
        let p_x = self.commitment.get().poly.clone();
        p_x - &self.low_degree_equivalent
    }
}

struct RotationSetExtension<'a, C: CurveAffine> {
    commitments: Vec<CommitmentExtension<'a, C>>,
    points: Vec<C::Scalar>,
}

impl<'a, C: CurveAffine> RotationSet<C::Scalar, PolynomialPointer<'a, C>> {
    fn extend(&self, commitments: Vec<CommitmentExtension<'a, C>>) -> RotationSetExtension<'a, C> {
        RotationSetExtension {
            commitments,
            points: self.points.clone(),
        }
    }
}

/// Concrete KZG prover with SHPLONK variant
#[derive(Debug)]
pub struct ProverSHPLONK<'a, E: Engine> {
    params: &'a ParamsKZG<E>,
}

impl<'a, E: Engine> ProverSHPLONK<'a, E> {
    /// Given parameters creates new prover instance
    pub fn new(params: &'a ParamsKZG<E>) -> Self {
        Self { params }
    }
}

/// Create a multi-opening proof
impl<'params, E: Engine + Debug> Prover<'params, KZGCommitmentScheme<E>>
    for ProverSHPLONK<'params, E>
{
    fn new(params: &'params ParamsKZG<E>) -> Self {
        Self { params }
    }

    /// Create a multi-opening proof
    fn create_proof<
        'com,
        Ch: EncodedChallenge<E::G1Affine>,
        T: TranscriptWrite<E::G1Affine, Ch>,
        R,
        I,
    >(
        &self,
        _: R,
        transcript: &mut T,
        queries: I,
    ) -> io::Result<()>
    where
        I: IntoIterator<Item = ProverQuery<'com, E::G1Affine>> + Clone,
        R: RngCore,
    {
        let zero = || Polynomial::<E::Scalar, Coeff> {
            values: vec![E::Scalar::zero(); self.params.n as usize],
            _marker: PhantomData,
        };

        // TODO: explore if it is safe to use same challenge
        // for different sets that are already combined with anoter challenge
        let y: ChallengeY<_> = transcript.squeeze_challenge_scalar();

        let quotient_contribution =
            |rotation_set: &RotationSetExtension<E::G1Affine>| -> Polynomial<E::Scalar, Coeff> {
                // [P_i_0(X) - R_i_0(X), P_i_1(X) - R_i_1(X), ... ]
                let numerators: Vec<Polynomial<E::Scalar, Coeff>> = rotation_set
                    .commitments
                    .iter()
                    .map(|commitment| commitment.quotient_contribution())
                    .collect();

                // define numerator polynomial as
                // N_i_j(X) = (P_i_j(X) - R_i_j(X))
                // and combine polynomials with same evaluation point set
                // N_i(X) = linear_combinination(y, N_i_j(X))
                // where y is random scalar to combine numerator polynomials
                let n_x: Polynomial<E::Scalar, Coeff> =
                    numerators.iter().fold(zero(), |acc, q_x| (acc * *y) + q_x);

                let points = &rotation_set.points[..];

                // quotient contribution of this evaluation set is
                // Q_i(X) = N_i(X) / Z_i(X) where
                // Z_i(X) = (x - r_i_0) * (x - r_i_1) * ...
                let mut poly = div_by_vanishing(n_x, points);
                poly.resize(self.params.n as usize, E::Scalar::zero());

                Polynomial {
                    values: poly,
                    _marker: PhantomData,
                }
            };

        let intermediate_sets = construct_intermediate_sets(queries);
        let (rotation_sets, super_point_set) = (
            intermediate_sets.rotation_sets,
            intermediate_sets.super_point_set,
        );

        let rotation_sets: Vec<RotationSetExtension<E::G1Affine>> = rotation_sets
            .iter()
            .map(|rotation_set| {
                let commitments: Vec<CommitmentExtension<E::G1Affine>> = rotation_set
                    .commitments
                    .iter()
                    .map(|commitment_data| {
                        commitment_data.extend(self.params.n, rotation_set.points.clone())
                    })
                    .collect();
                rotation_set.extend(commitments)
            })
            .collect();

        let v: ChallengeV<_> = transcript.squeeze_challenge_scalar();

        let quotient_polynomials: Vec<Polynomial<E::Scalar, Coeff>> =
            rotation_sets.iter().map(quotient_contribution).collect();

        let h_x: Polynomial<E::Scalar, Coeff> = quotient_polynomials
            .iter()
            .fold(zero(), |acc, u_x| (acc * *v) + u_x);

        let h = self.params.commit(&h_x, Blind::default()).to_affine();
        transcript.write_point(h)?;
        let u: ChallengeU<_> = transcript.squeeze_challenge_scalar();

        let zt_eval = evaluate_vanishing_polynomial(&super_point_set[..], *u);

        let linearisation_contribution =
            |rotation_set: RotationSetExtension<E::G1Affine>| -> (Polynomial<E::Scalar, Coeff>, E::Scalar) {
                let diffs: Vec<E::Scalar> = super_point_set
                    .iter()
                    .filter(|point| !rotation_set.points.contains(point))
                    .copied()
                    .collect();

                // calculate difference vanishing polynomial evaluation

                let z_i = evaluate_vanishing_polynomial(&diffs[..], *u);

                // inner linearisation contibutions are
                // [P_i_0(X) - r_i_0, P_i_1(X) - r_i_1, ... ] where
                // r_i_j = R_i_j(u) is the evaluation of low degree equivalent polynomial
                // where u is random evaluation point
                let inner_contributions: Vec<Polynomial<E::Scalar, Coeff>> = rotation_set
                    .commitments
                    .iter()
                    .map(|commitment| commitment.linearisation_contribution(*u))
                    .collect();

                // define inner contributor polynomial as
                // L_i_j(X) = (P_i_j(X) - r_i_j)
                // and combine polynomials with same evaluation point set
                // L_i(X) = linear_combinination(y, L_i_j(X))
                // where y is random scalar to combine inner contibutors
                let l_x: Polynomial<E::Scalar, Coeff> = inner_contributions
                    .iter()
                    .fold(zero(), |acc, l_x| (acc * *y) + l_x);

                // finally scale l_x by difference vanishing polynomial evaluation z_i
                (l_x * z_i, z_i)
            };

        #[allow(clippy::type_complexity)]
        let (linearisation_contibutions, z_diffs): (
            Vec<Polynomial<E::Scalar, Coeff>>,
            Vec<E::Scalar>,
        ) = rotation_sets
            .into_iter()
            .map(linearisation_contribution)
            .unzip();

        let l_x: Polynomial<E::Scalar, Coeff> = linearisation_contibutions
            .iter()
            .fold(zero(), |acc, u_x| (acc * *v) + u_x);

        let l_x = l_x - &(h_x * zt_eval);

        // sanity check
        {
            let must_be_zero = eval_polynomial(&l_x.values[..], *u);
            assert_eq!(must_be_zero, E::Scalar::zero());
        }

        let mut h_x = div_by_vanishing(l_x, &[*u]);

        // normalize coefficients by the coefficient of the first polynomial
        let z_0_diff_inv = z_diffs[0].invert().unwrap();
        for h_i in h_x.iter_mut() {
            h_i.mul_assign(z_0_diff_inv)
        }

        let h_x = Polynomial {
            values: h_x,
            _marker: PhantomData,
        };

        let h = self.params.commit(&h_x, Blind::default()).to_affine();
        transcript.write_point(h)?;

        Ok(())
    }
}
