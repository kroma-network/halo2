use std::marker::PhantomData;

use super::commitment::{IPACommitmentScheme, ParamsIPA, ParamsVerifierIPA};
use super::msm::MSMIPA;
use super::multiopen::VerifierIPA;
use crate::poly::commitment::CommitmentScheme;
use crate::transcript::TranscriptRead;
use crate::{
    arithmetic::best_multiexp,
    plonk::Error,
    poly::{
        commitment::MSM,
        strategy::{Guard, VerificationStrategy},
    },
    transcript::EncodedChallenge,
};
use ff::Field;
use group::Curve;
use halo2curves::CurveAffine;
use rand_core::RngCore;

/// Wrapper for verification accumulator
#[derive(Debug, Clone)]
pub struct GuardIPA<'params, C: CurveAffine> {
    pub(crate) msm_accumulator: MSMIPA<'params, C>,
    pub(crate) neg_c: C::Scalar,
    pub(crate) u: Vec<C::Scalar>,
    pub(crate) u_packed: Vec<C::Scalar>,
}

/// Define accumulator type as `MSMIPA`
impl<'params, C: CurveAffine> Guard<'params, IPACommitmentScheme<C>> for GuardIPA<'params, C> {
    type MSMAccumulator = MSMIPA<'params, C>;
}

/// IPA specific operations
impl<'params, C: CurveAffine> GuardIPA<'params, C> {
    /// Lets caller supply the challenges and obtain an MSM with updated
    /// scalars and points.
    pub fn use_challenges(mut self) -> MSMIPA<'params, C> {
        let s = compute_s(&self.u, self.neg_c);

        self.msm_accumulator.add_to_g_scalars(&s);

        self.msm_accumulator
    }

    /// Lets caller supply the purported G point and simply appends
    /// [-c] G to return an updated MSM.
    pub fn use_g(mut self, g: C) -> (MSMIPA<'params, C>, ProofAccumulator<C>) {
        self.msm_accumulator.append_term(self.neg_c, g.into());

        let accumulator = ProofAccumulator {
            g,
            u_packed: self.u_packed,
        };

        (self.msm_accumulator, accumulator)
    }

    /// Computes G = ⟨s, params.g⟩
    pub fn compute_g(&self) -> C {
        let s = compute_s(&self.u, C::Scalar::one());

        best_multiexp(&s, &self.msm_accumulator.params.g).to_affine()
    }
}

/// A verifier that checks multiple proofs in a batch.
#[derive(Debug)]
pub struct BatchVerifier<'params, C: CurveAffine, R: RngCore> {
    msm_accumulator: MSMIPA<'params, C>,
    rng: R,
}

impl<'params, C: CurveAffine, R: RngCore>
    VerificationStrategy<'params, IPACommitmentScheme<C>, VerifierIPA<'params, C>, R>
    for BatchVerifier<'params, C, R>
{
    type Output = Self;

    fn new(params: &'params ParamsIPA<C>, rng: R) -> Self {
        BatchVerifier {
            msm_accumulator: MSMIPA::new(params),
            rng,
        }
    }

    fn process(
        mut self,
        f: impl FnOnce(MSMIPA<'params, C>) -> Result<GuardIPA<'params, C>, Error>,
    ) -> Result<Self::Output, Error> {
        self.msm_accumulator.scale(C::Scalar::random(&mut self.rng));
        let guard = f(self.msm_accumulator)?;

        Ok(Self {
            msm_accumulator: guard.use_challenges(),
            rng: self.rng,
        })
    }

    /// Finalizes the batch and checks its validity.
    ///
    /// Returns `false` if *some* proof was invalid. If the caller needs to identify
    /// specific failing proofs, it must re-process the proofs separately.
    #[must_use]
    fn finalize(self) -> bool {
        self.msm_accumulator.check()
    }
}

/// Computes the coefficients of $g(X) = \prod\limits_{i=0}^{k-1} (1 + u_{k - 1 - i} X^{2^i})$.
fn compute_s<F: Field>(u: &[F], init: F) -> Vec<F> {
    assert!(!u.is_empty());
    let mut v = vec![F::zero(); 1 << u.len()];
    v[0] = init;

    for (len, u_j) in u.iter().rev().enumerate().map(|(i, u_j)| (1 << i, u_j)) {
        let (left, right) = v.split_at_mut(len);
        let right = &mut right[0..len];
        right.copy_from_slice(left);
        for v in right {
            *v *= u_j;
        }
    }

    v
}

/// An accumulator instance consisting of an evaluation claim and a proof.
#[derive(Debug, Clone)]
pub struct ProofAccumulator<C: CurveAffine> {
    /// The claimed output of the linear-time polycommit opening protocol
    pub g: C,

    /// A vector of challenges u_0, ..., u_{k - 1} sampled by the verifier, to
    /// be used in computing G'_0.
    pub u_packed: Vec<C::Scalar>,
}
