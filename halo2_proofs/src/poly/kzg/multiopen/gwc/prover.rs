use super::{construct_intermediate_sets, ChallengeV, Query};
use crate::arithmetic::{eval_polynomial, kate_division, CurveAffine, FieldExt};

use crate::poly::commitment::ParamsProver;
use crate::poly::commitment::Prover;
use crate::poly::kzg::commitment::{KZGCommitmentScheme, ParamsKZG};
use crate::poly::query::ProverQuery;
use crate::poly::Rotation;
use crate::poly::{
    commitment::{Blind, Params},
    Coeff, Polynomial,
};
use crate::transcript::{EncodedChallenge, TranscriptWrite};

use ff::Field;
use group::Curve;
use halo2curves::pairing::Engine;
use rand_core::RngCore;
use std::fmt::Debug;
use std::io::{self, Write};
use std::marker::PhantomData;

/// Concrete KZG prover with GWC variant
#[derive(Debug)]
pub struct ProverGWC<'params, E: Engine> {
    params: &'params ParamsKZG<E>,
}

/// Create a multi-opening proof
impl<'params, E: Engine + Debug, R: RngCore> Prover<'params, KZGCommitmentScheme<E>, R>
    for ProverGWC<'params, E>
{
    fn new(params: &'params ParamsKZG<E>) -> Self {
        Self { params }
    }

    /// Create a multi-opening proof
    fn create_proof<
        'com,
        Ch: EncodedChallenge<E::G1Affine>,
        T: TranscriptWrite<E::G1Affine, Ch>,
        I,
    >(
        &self,
        _: R,
        transcript: &mut T,
        queries: I,
    ) -> io::Result<()>
    where
        I: IntoIterator<Item = ProverQuery<'com, E::G1Affine>> + Clone,
    {
        let v: ChallengeV<_> = transcript.squeeze_challenge_scalar();
        let commitment_data = construct_intermediate_sets(queries);

        let zero = || Polynomial::<E::Scalar, Coeff> {
            values: vec![E::Scalar::zero(); self.params.n as usize],
            _marker: PhantomData,
        };

        for commitment_at_a_point in commitment_data.iter() {
            let mut poly_batch = zero();
            let mut eval_batch = E::Scalar::zero();
            let z = commitment_at_a_point.point;
            for query in commitment_at_a_point.queries.iter() {
                assert_eq!(query.get_point(), z);

                let poly = query.get_commitment().poly;
                let eval = query.get_eval();
                poly_batch = (poly_batch * *v) + poly;
                eval_batch = (eval_batch * *v) + eval;
            }

            let poly_batch = &poly_batch - eval_batch;
            let witness_poly = Polynomial {
                values: kate_division(&poly_batch.values, z),
                _marker: PhantomData,
            };
            let w = self
                .params
                .commit(&witness_poly, Blind::default())
                .to_affine();

            transcript.write_point(w)?;
        }
        Ok(())
    }
}
