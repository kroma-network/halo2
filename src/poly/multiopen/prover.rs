use super::super::{commitment::Params, Coeff, Polynomial};

use super::{construct_intermediate_sets, ChallengeV, ProverQuery, Query};

use crate::arithmetic::{eval_polynomial, kate_division, CurveAffine};
use crate::transcript::{EncodedChallenge, TranscriptWrite};

use ff::Field;
use group::Curve;
use pairing::arithmetic::Engine;
use std::io;
use std::marker::PhantomData;

#[derive(Debug, Clone)]
struct CommitmentData<C: CurveAffine> {
    set_index: usize,
    point_indices: Vec<usize>,
    evals: Vec<C::Scalar>,
}

/// Create a multi-opening proof
pub fn create_proof<
    'a,
    I,
    E: Engine,
    Ec: EncodedChallenge<E::G1Affine>,
    T: TranscriptWrite<E::G1Affine, Ec>,
>(
    params: &Params<E>,
    transcript: &mut T,
    queries: I,
) -> io::Result<()>
where
    I: IntoIterator<Item = ProverQuery<'a, E::G1Affine>> + Clone,
{
    let v: ChallengeV<_> = transcript.squeeze_challenge_scalar();
    let commitment_data = construct_intermediate_sets(queries);

    let zero = || Polynomial::<E::Scalar, Coeff> {
        values: vec![E::Scalar::zero(); params.n as usize],
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
            poly_batch = poly_batch * *v + poly;
            eval_batch = eval_batch * *v + eval;
        }

        let poly_batch = &poly_batch - eval_batch;
        let witness_poly = Polynomial {
            values: kate_division(&poly_batch.values, z),
            _marker: PhantomData,
        };
        let w = params.commit(&witness_poly).to_affine();
        transcript.write_point(w)?;
    }
    Ok(())
}

#[doc(hidden)]
#[derive(Copy, Clone)]
pub struct PolynomialPointer<'a, C: CurveAffine> {
    poly: &'a Polynomial<C::Scalar, Coeff>,
}

impl<'a, C: CurveAffine> PartialEq for PolynomialPointer<'a, C> {
    fn eq(&self, other: &Self) -> bool {
        std::ptr::eq(self.poly, other.poly)
    }
}

impl<'a, C: CurveAffine> Query<C::Scalar> for ProverQuery<'a, C> {
    type Commitment = PolynomialPointer<'a, C>;
    type Scalar = C::Scalar;

    fn get_point(&self) -> C::Scalar {
        self.point
    }
    fn get_eval(&self) -> C::Scalar {
        eval_polynomial(self.poly, self.get_point())
    }
    fn get_commitment(&self) -> Self::Commitment {
        PolynomialPointer { poly: self.poly }
    }
}
