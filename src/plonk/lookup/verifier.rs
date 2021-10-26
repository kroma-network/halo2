use std::iter;

use super::super::{
    circuit::Expression, ChallengeBeta, ChallengeGamma, ChallengeTheta, ChallengeX,
};
use super::Argument;
use crate::{
    arithmetic::{CurveAffine, FieldExt},
    plonk::{Error, VerifyingKey},
    poly::{multiopen::VerifierQuery, Rotation},
    transcript::{EncodedChallenge, TranscriptRead},
};
use ff::Field;
use pairing::arithmetic::Engine;

pub struct PermutationCommitments<E: Engine> {
    permuted_input_commitment: E::G1Affine,
    permuted_table_commitment: E::G1Affine,
}

pub struct Committed<E: Engine> {
    permuted: PermutationCommitments<E>,
    product_commitment: E::G1Affine,
}

pub struct Evaluated<E: Engine> {
    committed: Committed<E>,
    product_eval: E::Scalar,
    product_next_eval: E::Scalar,
    permuted_input_eval: E::Scalar,
    permuted_input_inv_eval: E::Scalar,
    permuted_table_eval: E::Scalar,
}

impl<F: FieldExt> Argument<F> {
    pub(in crate::plonk) fn read_permuted_commitments<
        E: Engine,
        Ec: EncodedChallenge<E::G1Affine>,
        T: TranscriptRead<E::G1Affine, Ec>,
    >(
        &self,
        transcript: &mut T,
    ) -> Result<PermutationCommitments<E>, Error> {
        let permuted_input_commitment = transcript
            .read_point()
            .map_err(|_| Error::TranscriptError)?;
        let permuted_table_commitment = transcript
            .read_point()
            .map_err(|_| Error::TranscriptError)?;

        Ok(PermutationCommitments {
            permuted_input_commitment,
            permuted_table_commitment,
        })
    }
}

impl<E: Engine> PermutationCommitments<E> {
    pub(in crate::plonk) fn read_product_commitment<
        Ec: EncodedChallenge<E::G1Affine>,
        T: TranscriptRead<E::G1Affine, Ec>,
    >(
        self,
        transcript: &mut T,
    ) -> Result<Committed<E>, Error> {
        let product_commitment = transcript
            .read_point()
            .map_err(|_| Error::TranscriptError)?;

        Ok(Committed {
            permuted: self,
            product_commitment,
        })
    }
}

impl<E: Engine> Committed<E> {
    pub(crate) fn evaluate<
        Ec: EncodedChallenge<E::G1Affine>,
        T: TranscriptRead<E::G1Affine, Ec>,
    >(
        self,
        transcript: &mut T,
    ) -> Result<Evaluated<E>, Error> {
        let product_eval = transcript
            .read_scalar()
            .map_err(|_| Error::TranscriptError)?;
        let product_next_eval = transcript
            .read_scalar()
            .map_err(|_| Error::TranscriptError)?;
        let permuted_input_eval = transcript
            .read_scalar()
            .map_err(|_| Error::TranscriptError)?;
        let permuted_input_inv_eval = transcript
            .read_scalar()
            .map_err(|_| Error::TranscriptError)?;
        let permuted_table_eval = transcript
            .read_scalar()
            .map_err(|_| Error::TranscriptError)?;

        Ok(Evaluated {
            committed: self,
            product_eval,
            product_next_eval,
            permuted_input_eval,
            permuted_input_inv_eval,
            permuted_table_eval,
        })
    }
}

impl<E: Engine> Evaluated<E> {
    pub(in crate::plonk) fn expressions<'a>(
        &'a self,
        l_0: E::Scalar,
        l_last: E::Scalar,
        l_blind: E::Scalar,
        argument: &'a Argument<E::Scalar>,
        theta: ChallengeTheta<E::G1Affine>,
        beta: ChallengeBeta<E::G1Affine>,
        gamma: ChallengeGamma<E::G1Affine>,
        advice_evals: &[E::Scalar],
        fixed_evals: &[E::Scalar],
        instance_evals: &[E::Scalar],
    ) -> impl Iterator<Item = E::Scalar> + 'a {
        let active_rows = E::Scalar::one() - (l_last + l_blind);

        let product_expression = || {
            // z(\omega X) (a'(X) + \beta) (s'(X) + \gamma)
            // - z(X) (\theta^{m-1} a_0(X) + ... + a_{m-1}(X) + \beta) (\theta^{m-1} s_0(X) + ... + s_{m-1}(X) + \gamma)
            let left = self.product_next_eval
                * &(self.permuted_input_eval + &*beta)
                * &(self.permuted_table_eval + &*gamma);

            let compress_expressions = |expressions: &[Expression<E::Scalar>]| {
                expressions
                    .iter()
                    .map(|expression| {
                        expression.evaluate(
                            &|scalar| scalar,
                            &|_| panic!("virtual selectors are removed during optimization"),
                            &|index, _, _| fixed_evals[index],
                            &|index, _, _| advice_evals[index],
                            &|index, _, _| instance_evals[index],
                            &|a| -a,
                            &|a, b| a + &b,
                            &|a, b| a * &b,
                            &|a, scalar| a * &scalar,
                        )
                    })
                    .fold(E::Scalar::zero(), |acc, eval| acc * &*theta + &eval)
            };
            let right = self.product_eval
                * &(compress_expressions(&argument.input_expressions) + &*beta)
                * &(compress_expressions(&argument.table_expressions) + &*gamma);

            (left - &right) * &active_rows
        };

        std::iter::empty()
            .chain(
                // l_0(X) * (1 - z'(X)) = 0
                Some(l_0 * &(E::Scalar::one() - &self.product_eval)),
            )
            .chain(
                // l_last(X) * (z(X)^2 - z(X)) = 0
                Some(l_last * &(self.product_eval.square() - &self.product_eval)),
            )
            .chain(
                // (1 - (l_last(X) + l_blind(X))) * (
                //   z(\omega X) (a'(X) + \beta) (s'(X) + \gamma)
                //   - z(X) (\theta^{m-1} a_0(X) + ... + a_{m-1}(X) + \beta) (\theta^{m-1} s_0(X) + ... + s_{m-1}(X) + \gamma)
                // ) = 0
                Some(product_expression()),
            )
            .chain(Some(
                // l_0(X) * (a'(X) - s'(X)) = 0
                l_0 * &(self.permuted_input_eval - &self.permuted_table_eval),
            ))
            .chain(Some(
                // (1 - (l_last(X) + l_blind(X))) * (a′(X) − s′(X))⋅(a′(X) − a′(\omega^{-1} X)) = 0
                (self.permuted_input_eval - &self.permuted_table_eval)
                    * &(self.permuted_input_eval - &self.permuted_input_inv_eval)
                    * &active_rows,
            ))
    }

    pub(in crate::plonk) fn queries<'r>(
        &'r self,
        vk: &'r VerifyingKey<E>,
        x: ChallengeX<E::G1Affine>,
    ) -> impl Iterator<Item = VerifierQuery<'r, E::G1Affine>> + Clone {
        let x_inv = vk.domain.rotate_omega(*x, Rotation::prev());
        let x_next = vk.domain.rotate_omega(*x, Rotation::next());

        iter::empty()
            // Open lookup product commitment at x
            .chain(Some(VerifierQuery::new_commitment(
                &self.committed.product_commitment,
                *x,
                self.product_eval,
            )))
            // Open lookup input commitments at x
            .chain(Some(VerifierQuery::new_commitment(
                &self.committed.permuted.permuted_input_commitment,
                *x,
                self.permuted_input_eval,
            )))
            // Open lookup table commitments at x
            .chain(Some(VerifierQuery::new_commitment(
                &self.committed.permuted.permuted_table_commitment,
                *x,
                self.permuted_table_eval,
            )))
            // Open lookup input commitments at \omega^{-1} x
            .chain(Some(VerifierQuery::new_commitment(
                &self.committed.permuted.permuted_input_commitment,
                x_inv,
                self.permuted_input_inv_eval,
            )))
            // Open lookup product commitment at \omega x
            .chain(Some(VerifierQuery::new_commitment(
                &self.committed.product_commitment,
                x_next,
                self.product_next_eval,
            )))
    }
}
