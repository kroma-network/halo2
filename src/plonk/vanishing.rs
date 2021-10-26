use std::marker::PhantomData;

use pairing::arithmetic::Engine;

use crate::arithmetic::CurveAffine;

mod prover;
mod verifier;

/// A vanishing argument.
pub(crate) struct Argument<E: Engine> {
    _marker: PhantomData<E>,
}
