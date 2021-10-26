//! This module provides an implementation of a variant of (Turbo)[PLONK][plonk]
//! that is designed specifically for the polynomial commitment scheme described
//! in the [Halo][halo] paper.
//!
//! [halo]: https://eprint.iacr.org/2019/1021
//! [plonk]: https://eprint.iacr.org/2019/953

use blake2b_simd::Params as Blake2bParams;
use group::GroupEncoding;
use pairing::arithmetic::Engine;

use crate::arithmetic::{BaseExt, CurveAffine, FieldExt};
use crate::poly::{
    commitment::Params, Coeff, EvaluationDomain, ExtendedLagrangeCoeff, LagrangeCoeff,
    PinnedEvaluationDomain, Polynomial,
};
use crate::transcript::{ChallengeScalar, EncodedChallenge, Transcript};

mod circuit;
mod keygen;
mod lookup;
pub(crate) mod permutation;
mod vanishing;

mod prover;
mod verifier;

pub use circuit::*;
pub use keygen::*;
pub use prover::*;
pub use verifier::*;

use std::io;

/// This is a verifying key which allows for the verification of proofs for a
/// particular circuit.
#[derive(Debug)]
pub struct VerifyingKey<E: Engine> {
    domain: EvaluationDomain<E::Scalar>,
    fixed_commitments: Vec<E::G1Affine>,
    permutation: permutation::VerifyingKey<E::G1Affine>,
    cs: ConstraintSystem<E::Scalar>,
}

impl<E: Engine> VerifyingKey<E> {
    /// Writes a verifying key to a buffer.
    pub fn write<W: io::Write>(&self, writer: &mut W) -> io::Result<()> {
        for commitment in &self.fixed_commitments {
            writer.write_all(commitment.to_bytes().as_ref())?;
        }
        self.permutation.write(writer)?;

        Ok(())
    }

    /// Reads a verification key from a buffer.
    pub fn read<R: io::Read, ConcreteCircuit: Circuit<E::Scalar>>(
        reader: &mut R,
        params: &Params<E>,
    ) -> io::Result<Self> {
        let (domain, cs, _) = keygen::create_domain::<E, ConcreteCircuit>(params);

        let fixed_commitments: Vec<_> = (0..cs.num_fixed_columns)
            .map(|_| E::G1Affine::read(reader))
            .collect::<Result<_, _>>()?;

        let permutation = permutation::VerifyingKey::read(reader, &cs.permutation)?;

        Ok(VerifyingKey {
            domain,
            fixed_commitments,
            permutation,
            cs,
        })
    }

    /// Hashes a verification key into a transcript.
    pub fn hash_into<Ec: EncodedChallenge<E::G1Affine>, T: Transcript<E::G1Affine, Ec>>(
        &self,
        transcript: &mut T,
    ) -> io::Result<()> {
        let mut hasher = Blake2bParams::new()
            .hash_length(64)
            .personal(b"Halo2-Verify-Key")
            .to_state();

        let s = format!("{:?}", self.pinned());

        hasher.update(&(s.len() as u64).to_le_bytes());
        hasher.update(s.as_bytes());

        // Hash in final Blake2bState
        transcript.common_scalar(E::Scalar::from_bytes_wide(hasher.finalize().as_array()))?;

        Ok(())
    }

    /// Obtains a pinned representation of this verification key that contains
    /// the minimal information necessary to reconstruct the verification key.
    pub fn pinned(&self) -> PinnedVerificationKey<'_, E::G1Affine> {
        PinnedVerificationKey {
            base_modulus: E::Scalar::MODULUS, // FIX
            scalar_modulus: E::Scalar::MODULUS,
            domain: self.domain.pinned(),
            fixed_commitments: &self.fixed_commitments,
            permutation: &self.permutation,
            cs: self.cs.pinned(),
        }
    }
}

/// Minimal representation of a verification key that can be used to identify
/// its active contents.
#[derive(Debug)]
pub struct PinnedVerificationKey<'a, C: CurveAffine> {
    base_modulus: &'static str,
    scalar_modulus: &'static str,
    domain: PinnedEvaluationDomain<'a, C::Scalar>,
    cs: PinnedConstraintSystem<'a, C::Scalar>,
    fixed_commitments: &'a Vec<C>,
    permutation: &'a permutation::VerifyingKey<C>,
}
/// This is a proving key which allows for the creation of proofs for a
/// particular circuit.
#[derive(Debug)]
pub struct ProvingKey<E: Engine> {
    vk: VerifyingKey<E>,
    l0: Polynomial<E::Scalar, ExtendedLagrangeCoeff>,
    l_blind: Polynomial<E::Scalar, ExtendedLagrangeCoeff>,
    l_last: Polynomial<E::Scalar, ExtendedLagrangeCoeff>,
    fixed_values: Vec<Polynomial<E::Scalar, LagrangeCoeff>>,
    fixed_polys: Vec<Polynomial<E::Scalar, Coeff>>,
    fixed_cosets: Vec<Polynomial<E::Scalar, ExtendedLagrangeCoeff>>,
    permutation: permutation::ProvingKey<E::G1Affine>,
}

/// This is an error that could occur during proving or circuit synthesis.
// TODO: these errors need to be cleaned up
#[derive(Debug, PartialEq)]
pub enum Error {
    /// This is an error that can occur during synthesis of the circuit, for
    /// example, when the witness is not present.
    SynthesisError,
    /// The structured reference string or the parameters are not compatible
    /// with the circuit being synthesized.
    IncompatibleParams,
    /// The constraint system is not satisfied.
    ConstraintSystemFailure,
    /// Out of bounds index passed to a backend
    BoundsFailure,
    /// Opening error
    OpeningError,
    /// Transcript error
    TranscriptError,
    /// Instance provided has more rows than supported by circuit
    NotEnoughRowsAvailable,
    /// Instance provided exceeds number of available rows
    InstanceTooLarge,
    /// Circuit synthesis requires global constants, but circuit configuration did not
    /// call [`ConstraintSystem::enable_constant`] on fixed columns with sufficient space.
    NotEnoughColumnsForConstants,
}

impl<E: Engine> ProvingKey<E> {
    /// Get the underlying [`VerifyingKey`].
    pub fn get_vk(&self) -> &VerifyingKey<E> {
        &self.vk
    }
}

impl<E: Engine> VerifyingKey<E> {
    /// Get the underlying [`EvaluationDomain`].
    pub fn get_domain(&self) -> &EvaluationDomain<E::Scalar> {
        &self.domain
    }
}

#[derive(Clone, Copy, Debug)]
struct Theta;
type ChallengeTheta<F> = ChallengeScalar<F, Theta>;

#[derive(Clone, Copy, Debug)]
struct Beta;
type ChallengeBeta<F> = ChallengeScalar<F, Beta>;

#[derive(Clone, Copy, Debug)]
struct Gamma;
type ChallengeGamma<F> = ChallengeScalar<F, Gamma>;

#[derive(Clone, Copy, Debug)]
struct Y;
type ChallengeY<F> = ChallengeScalar<F, Y>;

#[derive(Clone, Copy, Debug)]
struct X;
type ChallengeX<F> = ChallengeScalar<F, X>;
