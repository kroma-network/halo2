use group::Curve;
use pairing::arithmetic::Engine;
use std::iter;

use super::Argument;
use crate::{
    arithmetic::{eval_polynomial, BaseExt, CurveAffine},
    plonk::{ChallengeX, ChallengeY, Error},
    poly::{
        commitment::Params, multiopen::ProverQuery, Coeff, EvaluationDomain, ExtendedLagrangeCoeff,
        Polynomial,
    },
    transcript::{EncodedChallenge, TranscriptWrite},
};
use group::prime::PrimeCurveAffine;

pub(in crate::plonk) struct Committed<E: Engine> {
    random_poly: Polynomial<E::Scalar, Coeff>,
}

pub(in crate::plonk) struct Constructed<E: Engine> {
    h_pieces: Vec<Polynomial<E::Scalar, Coeff>>,
    committed: Committed<E>,
}

pub(in crate::plonk) struct Evaluated<E: Engine> {
    h_poly: Polynomial<E::Scalar, Coeff>,
    committed: Committed<E>,
}

impl<E: Engine> Argument<E> {
    pub(in crate::plonk) fn commit<
        Ec: EncodedChallenge<E::G1Affine>,
        T: TranscriptWrite<E::G1Affine, Ec>,
    >(
        params: &Params<E>,
        domain: &EvaluationDomain<E::Scalar>,
        transcript: &mut T,
    ) -> Result<Committed<E>, Error> {
        // Sample a random polynomial of degree n - 1
        let mut random_poly = domain.empty_coeff();
        for coeff in random_poly.iter_mut() {
            *coeff = E::Scalar::rand();
        }
        // Sample a random blinding factor

        // Commit
        let c = params.commit(&random_poly).to_affine();
        transcript
            .write_point(c)
            .map_err(|_| Error::TranscriptError)?;

        Ok(Committed { random_poly })
    }
}

impl<E: Engine> Committed<E> {
    pub(in crate::plonk) fn construct<
        Ec: EncodedChallenge<E::G1Affine>,
        T: TranscriptWrite<E::G1Affine, Ec>,
    >(
        self,
        params: &Params<E>,
        domain: &EvaluationDomain<E::Scalar>,
        expressions: impl Iterator<Item = Polynomial<E::Scalar, ExtendedLagrangeCoeff>>,
        y: ChallengeY<E::G1Affine>,
        transcript: &mut T,
    ) -> Result<Constructed<E>, Error> {
        // Evaluate the h(X) polynomial's constraint system expressions for the constraints provided
        let h_poly = expressions.fold(domain.empty_extended(), |h_poly, v| h_poly * *y + &v);

        // Divide by t(X) = X^{params.n} - 1.
        let h_poly = domain.divide_by_vanishing_poly(h_poly);

        // Obtain final h(X) polynomial
        let h_poly = domain.extended_to_coeff(h_poly);

        // Split h(X) up into pieces
        let h_pieces = h_poly
            .chunks_exact(params.n as usize)
            .map(|v| domain.coeff_from_vec(v.to_vec()))
            .collect::<Vec<_>>();
        drop(h_poly);

        // Compute commitments to each h(X) piece
        let h_commitments_projective: Vec<_> = h_pieces
            .iter()
            .map(|h_piece| params.commit(h_piece))
            .collect();
        let mut h_commitments = vec![E::G1Affine::identity(); h_commitments_projective.len()];
        E::G1::batch_normalize(&h_commitments_projective, &mut h_commitments);
        let h_commitments = h_commitments;

        // Hash each h(X) piece
        for c in h_commitments.iter() {
            transcript
                .write_point(*c)
                .map_err(|_| Error::TranscriptError)?;
        }

        Ok(Constructed {
            h_pieces,
            // h_blinds,
            committed: self,
        })
    }
}

impl<E: Engine> Constructed<E> {
    pub(in crate::plonk) fn evaluate<
        Ec: EncodedChallenge<E::G1Affine>,
        T: TranscriptWrite<E::G1Affine, Ec>,
    >(
        self,
        x: ChallengeX<E::G1Affine>,
        xn: E::Scalar,
        domain: &EvaluationDomain<E::Scalar>,
        transcript: &mut T,
    ) -> Result<Evaluated<E>, Error> {
        let h_poly = self
            .h_pieces
            .iter()
            .rev()
            .fold(domain.empty_coeff(), |acc, eval| acc * xn + eval);

        let random_eval = eval_polynomial(&self.committed.random_poly, *x);
        transcript
            .write_scalar(random_eval)
            .map_err(|_| Error::TranscriptError)?;

        Ok(Evaluated {
            h_poly,
            committed: self.committed,
        })
    }
}

impl<E: Engine> Evaluated<E> {
    pub(in crate::plonk) fn open(
        &self,
        x: ChallengeX<E::G1Affine>,
    ) -> impl Iterator<Item = ProverQuery<'_, E::G1Affine>> + Clone {
        iter::empty()
            .chain(Some(ProverQuery {
                point: *x,
                poly: &self.h_poly,
            }))
            .chain(Some(ProverQuery {
                point: *x,
                poly: &self.committed.random_poly,
            }))
    }
}
