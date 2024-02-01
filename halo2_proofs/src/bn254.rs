use std::{
    fmt,
    io::{self, Write},
    marker::PhantomData,
};

use crate::{
    plonk::{sealed, Column, Fixed},
    poly::commitment::Blind,
    transcript::{
        Challenge255, EncodedChallenge, Transcript, TranscriptWrite, TranscriptWriterBuffer,
    },
};
use ff::PrimeField;
use halo2curves::{Coordinates, CurveAffine};

#[repr(C)]
#[derive(Debug)]
pub struct G1Point2 {
    pub x: Fq,
    pub y: Fq,
}

#[repr(C)]
#[derive(Debug)]
pub struct G1JacobianPoint {
    pub x: Fq,
    pub y: Fq,
    pub z: Fq,
}

#[repr(transparent)]
#[derive(Debug)]
pub struct Fq(pub [u64; 4]);

#[repr(transparent)]
#[derive(Debug)]
pub struct Fr(pub [u64; 4]);

#[derive(Debug)]
pub struct InstanceSingle {
    pub instance_values: Vec<Evals>,
    pub instance_polys: Vec<Poly>,
}

#[derive(Clone, Debug)]
pub struct AdviceSingle {
    pub advice_polys: Vec<Evals>,
    pub advice_blinds: Vec<Blind<halo2curves::bn256::Fr>>,
}

#[cxx::bridge(namespace = "tachyon::halo2_api::bn254")]
pub mod ffi {
    extern "Rust" {
        type G1JacobianPoint;
        type G1Point2;
        type Fr;
        type InstanceSingle;
        type AdviceSingle;
    }

    unsafe extern "C++" {
        include!("halo2_proofs/include/bn254_blake2b_writer.h");

        type Blake2bWriter;

        fn new_blake2b_writer() -> UniquePtr<Blake2bWriter>;
        fn update(self: Pin<&mut Blake2bWriter>, data: &[u8]);
        fn finalize(self: Pin<&mut Blake2bWriter>, result: &mut [u8; 64]);
        fn state(&self) -> Vec<u8>;
    }

    unsafe extern "C++" {
        include!("halo2_proofs/include/bn254_proving_key.h");

        type ProvingKey;

        fn new_proving_key(data: &[u8]) -> UniquePtr<ProvingKey>;
        fn advice_column_phases(&self) -> Vec<u8>;
        fn blinding_factors(&self) -> u32;
        fn challenge_phases(&self) -> Vec<u8>;
        fn constants(&self) -> Vec<usize>;
        fn num_advice_columns(&self) -> usize;
        fn num_challenges(&self) -> usize;
        fn num_instance_columns(&self) -> usize;
        fn phases(&self) -> Vec<u8>;
        fn transcript_repr(self: Pin<&mut ProvingKey>, prover: &SHPlonkProver) -> Box<Fr>;
    }

    unsafe extern "C++" {
        include!("halo2_proofs/include/bn254_evals.h");

        type Evals;

        fn zero_evals() -> UniquePtr<Evals>;
        fn len(&self) -> usize;
        fn set_value(self: Pin<&mut Evals>, idx: usize, value: &Fr);
        fn clone(&self) -> UniquePtr<Evals>;
    }

    unsafe extern "C++" {
        include!("halo2_proofs/include/bn254_rational_evals.h");

        type RationalEvals;

        fn set_zero(self: Pin<&mut RationalEvals>, idx: usize);
        fn set_trivial(self: Pin<&mut RationalEvals>, idx: usize, numerator: &Fr);
        fn set_rational(
            self: Pin<&mut RationalEvals>,
            idx: usize,
            numerator: &Fr,
            denominator: &Fr,
        );
        fn clone(&self) -> UniquePtr<RationalEvals>;
    }

    unsafe extern "C++" {
        include!("halo2_proofs/include/bn254_poly.h");

        type Poly;
    }

    unsafe extern "C++" {
        include!("halo2_proofs/include/bn254_shplonk_prover.h");

        type SHPlonkProver;

        fn new_shplonk_prover(k: u32, s: &Fr) -> UniquePtr<SHPlonkProver>;
        fn k(&self) -> u32;
        fn n(&self) -> u64;
        fn commit(&self, poly: &Poly) -> Box<G1JacobianPoint>;
        fn commit_lagrange(&self, evals: &Evals) -> Box<G1JacobianPoint>;
        fn empty_evals(&self) -> UniquePtr<Evals>;
        fn empty_rational_evals(&self) -> UniquePtr<RationalEvals>;
        fn ifft(&self, evals: &Evals) -> UniquePtr<Poly>;
        fn batch_evaluate(
            &self,
            rational_evals: &[UniquePtr<RationalEvals>],
            evals: &mut [UniquePtr<Evals>],
        );
        fn set_rng(self: Pin<&mut SHPlonkProver>, state: &[u8]);
        fn set_transcript(self: Pin<&mut SHPlonkProver>, state: &[u8]);
        fn set_extended_domain(self: Pin<&mut SHPlonkProver>, pk: &ProvingKey);
        fn create_proof(
            self: Pin<&mut SHPlonkProver>,
            key: &ProvingKey,
            instance_singles: &mut [InstanceSingle],
            advice_singles: &mut [AdviceSingle],
            challenges: &[Fr],
        );
        fn get_proof(self: &SHPlonkProver) -> Vec<u8>;
    }
}

impl fmt::Debug for ffi::Blake2bWriter {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Blake2bWriter").finish()
    }
}

impl fmt::Debug for ffi::ProvingKey {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("ProvingKey").finish()
    }
}

impl fmt::Debug for ffi::Evals {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Evals").finish()
    }
}

impl fmt::Debug for ffi::RationalEvals {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("RationalEvals").finish()
    }
}

impl fmt::Debug for ffi::Poly {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Poly").finish()
    }
}

impl fmt::Debug for ffi::SHPlonkProver {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("SHPlonkProver").finish()
    }
}

#[derive(Debug)]
pub struct Blake2bWrite<W: Write, C: CurveAffine, E: EncodedChallenge<C>> {
    state: cxx::UniquePtr<ffi::Blake2bWriter>,
    writer: W,
    _marker: PhantomData<(W, C, E)>,
}

impl<W: Write, C: CurveAffine> Blake2bWrite<W, C, Challenge255<C>> {
    pub fn state(&self) -> Vec<u8> {
        self.state.state()
    }
}

impl<W: Write, C: CurveAffine> Transcript<C, Challenge255<C>>
    for Blake2bWrite<W, C, Challenge255<C>>
{
    fn squeeze_challenge(&mut self) -> Challenge255<C> {
        // Prefix to a prover's message soliciting a challenge
        const BLAKE2B_PREFIX_CHALLENGE: u8 = 0;
        self.state.pin_mut().update(&[BLAKE2B_PREFIX_CHALLENGE]);
        let mut result: [u8; 64] = [0; 64];
        self.state.pin_mut().finalize(&mut result);
        Challenge255::<C>::new(&result)
    }

    fn common_point(&mut self, point: C) -> io::Result<()> {
        // Prefix to a prover's message containing a curve point
        const BLAKE2B_PREFIX_POINT: u8 = 1;
        self.state.pin_mut().update(&[BLAKE2B_PREFIX_POINT]);
        let coords: Coordinates<C> = Option::from(point.coordinates()).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::Other,
                "cannot write points at infinity to the transcript",
            )
        })?;
        self.state.pin_mut().update(coords.x().to_repr().as_ref());
        self.state.pin_mut().update(coords.y().to_repr().as_ref());

        Ok(())
    }

    fn common_scalar(&mut self, scalar: C::Scalar) -> io::Result<()> {
        // Prefix to a prover's message containing a scalar
        const BLAKE2B_PREFIX_SCALAR: u8 = 2;
        self.state.pin_mut().update(&[BLAKE2B_PREFIX_SCALAR]);
        self.state.pin_mut().update(scalar.to_repr().as_ref());
        Ok(())
    }
}

impl<W: Write, C: CurveAffine> TranscriptWrite<C, Challenge255<C>>
    for Blake2bWrite<W, C, Challenge255<C>>
{
    fn write_point(&mut self, point: C) -> io::Result<()> {
        self.common_point(point)?;
        let compressed = point.to_bytes();
        self.writer.write_all(compressed.as_ref())
    }

    fn write_scalar(&mut self, scalar: C::Scalar) -> io::Result<()> {
        self.common_scalar(scalar)?;
        let data = scalar.to_repr();
        self.writer.write_all(data.as_ref())
    }
}

impl<W: Write, C: CurveAffine> TranscriptWriterBuffer<W, C, Challenge255<C>>
    for Blake2bWrite<W, C, Challenge255<C>>
{
    /// Initialize a transcript given an output buffer.
    fn init(writer: W) -> Self {
        Blake2bWrite {
            state: ffi::new_blake2b_writer(),
            writer: writer,
            _marker: PhantomData,
        }
    }

    fn finalize(self) -> W {
        // TODO: handle outstanding scalars? see issue #138
        self.writer
    }
}

#[derive(Debug)]
pub struct ProvingKey {
    inner: cxx::UniquePtr<ffi::ProvingKey>,
}

impl ProvingKey {
    pub fn from(data: &[u8]) -> ProvingKey {
        ProvingKey {
            inner: ffi::new_proving_key(data),
        }
    }

    // NOTE(chokobole): We name this as plural since it contains multi phases.
    // pk.vk.cs.advice_column_phase
    pub fn advice_column_phases(&self) -> Vec<sealed::Phase> {
        unsafe {
            let phases: Vec<sealed::Phase> = std::mem::transmute(self.inner.advice_column_phases());
            phases
        }
    }

    // pk.vk.cs.blinding_factors()
    pub fn blinding_factors(&self) -> u32 {
        self.inner.blinding_factors()
    }

    // NOTE(chokobole): We name this as plural since it contains multi phases.
    // pk.vk.cs.challenge_phase
    pub fn challenge_phases(&self) -> Vec<sealed::Phase> {
        unsafe {
            let phases: Vec<sealed::Phase> = std::mem::transmute(self.inner.challenge_phases());
            phases
        }
    }

    // pk.vk.cs.constants
    pub fn constants(&self) -> Vec<Column<Fixed>> {
        let constants = self
            .inner
            .constants()
            .iter()
            .map(|index| Column {
                index: *index,
                column_type: Fixed,
            })
            .collect::<Vec<_>>();
        constants
    }

    // pk.vk.cs.num_advice_columns
    pub fn num_advice_columns(&self) -> usize {
        self.inner.num_advice_columns()
    }

    // pk.vk.cs.num_challenges
    pub fn num_challenges(&self) -> usize {
        self.inner.num_challenges()
    }

    // pk.vk.cs.num_instance_columns
    pub fn num_instance_columns(&self) -> usize {
        self.inner.num_instance_columns()
    }

    // pk.vk.cs.phases()
    pub fn phases(&self) -> Vec<sealed::Phase> {
        unsafe {
            let phases: Vec<sealed::Phase> = std::mem::transmute(self.inner.phases());
            phases
        }
    }

    // pk.vk.transcript_repr
    pub fn transcript_repr(&mut self, prover: &SHPlonkProver) -> halo2curves::bn256::Fr {
        *unsafe {
            std::mem::transmute::<_, Box<halo2curves::bn256::Fr>>(
                self.inner.pin_mut().transcript_repr(&prover.inner),
            )
        }
    }
}

#[derive(Debug)]
pub struct Evals {
    inner: cxx::UniquePtr<ffi::Evals>,
}

impl Evals {
    pub fn zero() -> Evals {
        Self::new(ffi::zero_evals())
    }

    pub fn new(inner: cxx::UniquePtr<ffi::Evals>) -> Evals {
        Evals { inner }
    }

    pub fn len(&self) -> usize {
        self.inner.len()
    }

    pub fn set_value(&mut self, idx: usize, fr: &halo2curves::bn256::Fr) {
        let cpp_fr = unsafe { std::mem::transmute::<_, &Fr>(fr) };
        self.inner.pin_mut().set_value(idx, cpp_fr)
    }
}

impl Clone for Evals {
    fn clone(&self) -> Self {
        Self {
            inner: self.inner.clone(),
        }
    }
}

#[derive(Debug)]
pub struct RationalEvals {
    inner: cxx::UniquePtr<ffi::RationalEvals>,
}

impl RationalEvals {
    pub fn new(inner: cxx::UniquePtr<ffi::RationalEvals>) -> RationalEvals {
        RationalEvals { inner }
    }

    pub fn set_zero(&mut self, idx: usize) {
        self.inner.pin_mut().set_zero(idx)
    }

    pub fn set_trivial(&mut self, idx: usize, numerator: &halo2curves::bn256::Fr) {
        let cpp_numerator = unsafe { std::mem::transmute::<_, &Fr>(numerator) };
        self.inner.pin_mut().set_trivial(idx, cpp_numerator)
    }

    pub fn set_rational(
        &mut self,
        idx: usize,
        numerator: &halo2curves::bn256::Fr,
        denominator: &halo2curves::bn256::Fr,
    ) {
        let cpp_numerator = unsafe { std::mem::transmute::<_, &Fr>(numerator) };
        let cpp_denominator = unsafe { std::mem::transmute::<_, &Fr>(denominator) };
        self.inner
            .pin_mut()
            .set_rational(idx, cpp_numerator, cpp_denominator)
    }
}

impl Clone for RationalEvals {
    fn clone(&self) -> Self {
        Self {
            inner: self.inner.clone(),
        }
    }
}

#[derive(Debug)]
pub struct Poly {
    inner: cxx::UniquePtr<ffi::Poly>,
}

impl Poly {
    pub fn new(inner: cxx::UniquePtr<ffi::Poly>) -> Poly {
        Poly { inner }
    }
}

#[derive(Debug)]
pub struct SHPlonkProver {
    inner: cxx::UniquePtr<ffi::SHPlonkProver>,
}

impl SHPlonkProver {
    pub fn new(k: u32, s: &halo2curves::bn256::Fr) -> SHPlonkProver {
        let cpp_s = unsafe { std::mem::transmute::<_, &Fr>(s) };
        SHPlonkProver {
            inner: ffi::new_shplonk_prover(k, cpp_s),
        }
    }

    pub fn k(&self) -> u32 {
        self.inner.k()
    }

    pub fn n(&self) -> u64 {
        self.inner.n()
    }

    pub fn commit(&self, poly: &Poly) -> halo2curves::bn256::G1 {
        *unsafe {
            std::mem::transmute::<_, Box<halo2curves::bn256::G1>>(self.inner.commit(&poly.inner))
        }
    }

    pub fn commit_lagrange(&self, evals: &Evals) -> halo2curves::bn256::G1 {
        *unsafe {
            std::mem::transmute::<_, Box<halo2curves::bn256::G1>>(
                self.inner.commit_lagrange(&evals.inner),
            )
        }
    }

    pub fn empty_evals(&self) -> Evals {
        Evals::new(self.inner.empty_evals())
    }

    pub fn empty_rational_evals(&self) -> RationalEvals {
        RationalEvals::new(self.inner.empty_rational_evals())
    }

    pub fn batch_evaluate(&self, rational_evals: &[RationalEvals], evals: &mut [Evals]) {
        unsafe {
            let rational_evals: &[cxx::UniquePtr<ffi::RationalEvals>] =
                std::mem::transmute(rational_evals);
            let evals: &mut [cxx::UniquePtr<ffi::Evals>] = std::mem::transmute(evals);
            self.inner.batch_evaluate(rational_evals, evals)
        }
    }

    pub fn ifft(&self, evals: &Evals) -> Poly {
        Poly::new(self.inner.ifft(&evals.inner))
    }

    pub fn set_rng(&mut self, state: &[u8]) {
        self.inner.pin_mut().set_rng(state)
    }

    pub fn set_transcript(&mut self, state: &[u8]) {
        self.inner.pin_mut().set_transcript(state)
    }

    pub fn set_extended_domain(&mut self, pk: &ProvingKey) {
        self.inner.pin_mut().set_extended_domain(&pk.inner)
    }

    pub fn create_proof(
        &mut self,
        key: &ProvingKey,
        instance_singles: &mut [InstanceSingle],
        advice_singles: &mut [AdviceSingle],
        challenges: &[Fr],
    ) {
        self.inner
            .pin_mut()
            .create_proof(&key.inner, instance_singles, advice_singles, challenges)
    }

    pub fn get_proof(&self) -> Vec<u8> {
        self.inner.get_proof()
    }
}
