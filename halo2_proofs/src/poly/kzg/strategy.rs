use std::{fmt::Debug, marker::PhantomData};

use super::{
    commitment::{KZGCommitmentScheme, ParamsKZG},
    msm::{DualMSM, MSMKZG},
    multiopen::VerifierGWC,
};
use crate::{
    plonk::Error,
    poly::{
        commitment::{Verifier, MSM},
        ipa::msm::MSMIPA,
        strategy::{Guard, VerificationStrategy},
    },
    transcript::{EncodedChallenge, TranscriptRead},
};
use ff::Field;
use group::Group;
use halo2curves::{
    pairing::{Engine, MillerLoopResult, MultiMillerLoop},
    CurveAffine,
};
use rand_core::RngCore;

/// Wrapper for linear verification accumulator
#[derive(Debug, Clone)]
pub struct GuardKZG<'params, E: MultiMillerLoop + Debug> {
    pub(crate) msm_accumulator: DualMSM<'params, E>,
}

/// Define accumulator type as `DualMSM`
impl<'params, E: MultiMillerLoop + Debug> Guard<'params, KZGCommitmentScheme<E>>
    for GuardKZG<'params, E>
{
    type MSMAccumulator = DualMSM<'params, E>;
}

/// KZG specific operations
impl<'params, E: MultiMillerLoop + Debug> GuardKZG<'params, E> {
    pub(crate) fn new(msm_accumulator: DualMSM<'params, E>) -> Self {
        Self { msm_accumulator }
    }
}

/// A verifier that checks multiple proofs in a batch
#[derive(Clone, Debug)]
pub struct BatchVerifier<'params, E: Engine, R: RngCore> {
    pub(crate) msm_accumulator: DualMSM<'params, E>,
    pub(crate) rng: R,
}

impl<'params, E: MultiMillerLoop + Debug, R: RngCore> BatchVerifier<'params, E, R> {
    /// Constructs an empty batch verifier
    pub fn new(params: &'params ParamsKZG<E>, rng: R) -> Self {
        BatchVerifier {
            msm_accumulator: DualMSM::new(params),
            rng,
        }
    }

    /// Constructs and initialized new batch verifier
    pub fn with(msm_accumulator: DualMSM<'params, E>, rng: R) -> Self {
        BatchVerifier {
            msm_accumulator,
            rng,
        }
    }
}
