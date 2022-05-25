use std::fmt::Debug;
use std::io::Read;

use super::ChallengeY;
use super::{construct_intermediate_sets, ChallengeU, ChallengeV};
use crate::arithmetic::{
    eval_polynomial, evaluate_vanishing_polynomial, lagrange_interpolate, CurveAffine, FieldExt,
};
use crate::poly::commitment::Verifier;
use crate::poly::commitment::MSM;
use crate::poly::kzg::commitment::{KZGCommitmentScheme, ParamsKZG};
use crate::poly::kzg::msm::DualMSM;
use crate::poly::kzg::msm::{PreMSM, MSMKZG};
use crate::poly::kzg::strategy::{BatchVerifier, GuardKZG};
use crate::poly::query::Query;
use crate::poly::query::{CommitmentReference, VerifierQuery};
use crate::poly::strategy::VerificationStrategy;
use crate::poly::{
    commitment::{Params, ParamsVerifier},
    Error,
};
use crate::transcript::{EncodedChallenge, TranscriptRead};
use ff::Field;
use group::Group;
use halo2curves::pairing::{Engine, MillerLoopResult, MultiMillerLoop};
use rand_core::RngCore;
use std::ops::MulAssign;

/// Concrete KZG multiopen verifier with SHPLONK variant
#[derive(Debug)]
pub struct VerifierSHPLONK<'params, E: Engine> {
    params: &'params ParamsKZG<E>,
}

impl<'params, E: MultiMillerLoop + Debug> Verifier<'params, KZGCommitmentScheme<E>>
    for VerifierSHPLONK<'params, E>
{
    type Guard = GuardKZG<'params, E>;
    type MSMAccumulator = DualMSM<'params, E>;

    fn new(params: &'params ParamsKZG<E>) -> Self {
        Self { params }
    }

    /// Verify a multi-opening proof
    fn verify_proof<
        'com,
        Ch: EncodedChallenge<E::G1Affine>,
        T: TranscriptRead<E::G1Affine, Ch>,
        I,
    >(
        &self,
        transcript: &mut T,
        queries: I,
        mut msm_accumulator: DualMSM<'params, E>,
    ) -> Result<Self::Guard, Error>
    where
        I: IntoIterator<Item = VerifierQuery<'com, E::G1Affine>> + Clone,
    {
        let intermediate_sets = construct_intermediate_sets(queries);
        let (rotation_sets, super_point_set) = (
            intermediate_sets.rotation_sets,
            intermediate_sets.super_point_set,
        );

        let y: ChallengeY<_> = transcript.squeeze_challenge_scalar();
        let v: ChallengeV<_> = transcript.squeeze_challenge_scalar();

        let h1 = transcript.read_point().map_err(|_| Error::SamplingError)?;
        let u: ChallengeU<_> = transcript.squeeze_challenge_scalar();
        let h2 = transcript.read_point().map_err(|_| Error::SamplingError)?;

        let (mut z_0_diff_inverse, mut z_0) = (E::Scalar::zero(), E::Scalar::zero());
        let (mut outer_msm, mut r_outer_acc) = (PreMSM::<E>::new(), E::Scalar::zero());
        for (i, rotation_set) in rotation_sets.iter().enumerate() {
            let diffs: Vec<E::Scalar> = super_point_set
                .iter()
                .filter(|point| !rotation_set.points.contains(point))
                .copied()
                .collect();
            let mut z_diff_i = evaluate_vanishing_polynomial(&diffs[..], *u);

            // normalize coefficients by the coefficient of the first commitment
            if i == 0 {
                z_0 = evaluate_vanishing_polynomial(&rotation_set.points[..], *u);
                z_0_diff_inverse = z_diff_i.invert().unwrap();
                z_diff_i = E::Scalar::one();
            } else {
                z_diff_i.mul_assign(z_0_diff_inverse);
            }

            let (mut inner_msm, mut r_inner_acc) = (MSMKZG::new(), E::Scalar::zero());
            for commitment_data in rotation_set.commitments.iter() {
                // calculate low degree equivalent
                let r_x =
                    lagrange_interpolate(&rotation_set.points[..], &commitment_data.evals()[..]);
                let r_eval = eval_polynomial(&r_x[..], *u);
                r_inner_acc = (*y * r_inner_acc) + r_eval;

                let inner_contrib = match commitment_data.get() {
                    CommitmentReference::Commitment(c) => (*c).into(),
                    // TODO: we should support one more nested degree to append
                    // folded commitments to the inner_msm
                    CommitmentReference::MSM(msm) => msm.eval(),
                };
                inner_msm.append_term(E::Scalar::one(), inner_contrib);
            }
            r_outer_acc = (*v * r_outer_acc) + (r_inner_acc * z_diff_i);

            inner_msm.combine_with_base(*y);
            inner_msm.scale(z_diff_i);
            outer_msm.add_msm(inner_msm);
        }
        outer_msm.combine_with_base(*v);
        let mut outer_msm = outer_msm.normalize();
        let g1: E::G1 = self.params.g[0].into();
        outer_msm.append_term(-r_outer_acc, g1);
        outer_msm.append_term(-z_0, h1.into());
        outer_msm.append_term(*u, h2.into());

        msm_accumulator
            .left
            .append_term(E::Scalar::one(), h2.into());

        msm_accumulator.right.add_msm(&outer_msm);

        Ok(Self::Guard::new(msm_accumulator))
    }
}

impl<'params, E: MultiMillerLoop + Debug, R: RngCore>
    VerificationStrategy<'params, KZGCommitmentScheme<E>, VerifierSHPLONK<'params, E>, R>
    for BatchVerifier<'params, E, R>
{
    type Output = Self;

    /// Constructs a new batch verifier.
    fn new(params: &'params ParamsKZG<E>, rng: R) -> Self {
        BatchVerifier::new(params, rng)
    }

    fn process(
        mut self,

        f: impl FnOnce(DualMSM<'params, E>) -> Result<GuardKZG<'params, E>, crate::plonk::Error>,
    ) -> Result<Self::Output, crate::plonk::Error> {
        self.msm_accumulator.scale(E::Scalar::random(&mut self.rng));

        // Guard is updated with new msm contributions
        let guard = f(self.msm_accumulator)?;
        Ok(BatchVerifier::with(guard.msm_accumulator, self.rng))
    }

    fn finalize(self) -> bool {
        self.msm_accumulator.check()
    }
}
