//! This module contains an implementation of the polynomial commitment scheme
//! described in the [Halo][halo] paper.
//!
//! [halo]: https://eprint.iacr.org/2019/1021

use super::{Coeff, FFTData, LagrangeCoeff, Polynomial, MSM};
use crate::arithmetic::{
    best_fft, best_multiexp, parallelize, CurveAffine, CurveExt, Engine, FieldExt, Group,
};
use crate::helpers::CurveRead;
use crate::poly::FFTType;

use ff::{Field, PrimeField};
use group::{prime::PrimeCurveAffine, Curve, Group as _, GroupEncoding};
use rand_core::OsRng;
use rayon::{current_num_threads, scope};
use std::io;
use std::marker::PhantomData;
use std::ops::{Add, AddAssign, Mul, MulAssign};

/// These are the prover parameters for the polynomial commitment scheme.
#[derive(Debug)]
pub struct Params<C: CurveAffine> {
    pub(crate) k: u32,
    pub(crate) n: u64,
    pub(crate) g: Vec<C>,
    pub(crate) g_lagrange: Vec<C>,
    pub(crate) additional_data: Vec<u8>,
}

/// These are the verifier parameters for the polynomial commitment scheme.
#[derive(Debug)]
pub struct ParamsVerifier<E: Engine> {
    pub(crate) k: u32,
    pub(crate) n: u64,
    pub(crate) g1: E::G1Affine,
    pub(crate) g2: E::G2Affine,
    pub(crate) s_g2: E::G2Affine,
    pub(crate) g_lagrange: Vec<E::G1Affine>,
}

#[cfg(test)]
mod tests {
    use pairing::arithmetic::CurveAffine;

    use super::Params;

    impl<C: CurveAffine> Params<C> {
        pub fn do_this() {}
    }
}

impl<C: CurveAffine> Params<C> {
    /// Initializes parameters for the curve, Draws random toxic point inside of the function
    /// MUST NOT be used in production
    pub fn unsafe_setup<E: Engine>(k: u32) -> Params<E::G1Affine> {
        // TODO: Make this function only available in test mod
        // Largest root of unity exponent of the Engine is `2^E::Scalar::S`, so we can
        // only support FFTs of polynomials below degree `2^E::Scalar::S`.
        assert!(k <= E::Scalar::S);
        let n: u64 = 1 << k;

        let s = E::Scalar::random(OsRng);

        let mut g_projective: Vec<E::G1> = Vec::with_capacity(n as usize);
        let g1 = <E::G1Affine as PrimeCurveAffine>::generator();
        g_projective.push(g1.into());
        // g = [G1, [s] G1, [s^2] G1, ..., [s^(n-1)] G1]
        for i in 1..(n as usize) {
            g_projective.push(g_projective[i - 1] * s);
        }

        let g = {
            let mut g = vec![E::G1Affine::identity(); n as usize];
            parallelize(&mut g, |g, starts| {
                E::G1::batch_normalize(&g_projective[starts..(starts + g.len())], g);
            });
            g
        };

        // Let's evaluate all of the Lagrange basis polynomials
        // using an inverse FFT.
        let mut alpha_inv = E::Scalar::ROOT_OF_UNITY_INV;
        for _ in k..E::Scalar::S {
            alpha_inv = alpha_inv.square();
        }
        let mut g_lagrange_projective = g_projective;
        prev_fft(&mut g_lagrange_projective, alpha_inv, k);
        let minv = E::Scalar::TWO_INV.pow_vartime(&[k as u64, 0, 0, 0]);
        parallelize(&mut g_lagrange_projective, |g, _| {
            for g in g.iter_mut() {
                *g *= minv;
            }
        });

        let g_lagrange = {
            let mut g_lagrange = vec![E::G1Affine::identity(); n as usize];
            parallelize(&mut g_lagrange, |g_lagrange, starts| {
                E::G1::batch_normalize(
                    &g_lagrange_projective[starts..(starts + g_lagrange.len())],
                    g_lagrange,
                );
            });
            drop(g_lagrange_projective);
            g_lagrange
        };

        let g2 = <E::G2Affine as PrimeCurveAffine>::generator();
        let s_g2 = g2 * s;
        let additional_data = Vec::from(s_g2.to_bytes().as_ref());
        Params {
            k,
            n,
            g,
            g_lagrange,
            additional_data,
        }
    }

    /// This computes a commitment to a polynomial described by the provided
    /// slice of coefficients. The commitment will be blinded by the blinding
    /// factor `r`.
    pub fn commit(&self, poly: &Polynomial<C::Scalar, Coeff>) -> C::Curve {
        let mut scalars = Vec::with_capacity(poly.len());
        scalars.extend(poly.iter());
        let bases = &self.g;
        let size = scalars.len();
        assert!(bases.len() >= size);
        best_multiexp(&scalars, &bases[0..size])
    }

    /// This commits to a polynomial using its evaluations over the $2^k$ size
    /// evaluation domain. The commitment will be blinded by the blinding factor
    /// `r`.
    pub fn commit_lagrange(&self, poly: &Polynomial<C::Scalar, LagrangeCoeff>) -> C::Curve {
        let mut scalars = Vec::with_capacity(poly.len());
        scalars.extend(poly.iter());
        let bases = &self.g_lagrange;
        let size = scalars.len();
        assert!(bases.len() >= size);
        best_multiexp(&scalars, &bases[0..size])
    }

    /// Generates an empty multiscalar multiplication struct using the
    /// appropriate params.
    pub fn empty_msm(&self) -> MSM<C> {
        MSM::new()
    }

    /// Getter for g generators
    pub fn get_g(&self) -> Vec<C> {
        self.g.clone()
    }

    /// Writes params to a buffer.
    pub fn write<W: io::Write>(&self, writer: &mut W) -> io::Result<()> {
        writer.write_all(&self.k.to_le_bytes())?;
        for el in &self.g {
            writer.write_all(el.to_bytes().as_ref())?;
        }
        for el in &self.g_lagrange {
            writer.write_all(el.to_bytes().as_ref())?;
        }
        let additional_data_len = self.additional_data.len() as u32;
        writer.write_all(&additional_data_len.to_le_bytes())?;
        writer.write_all(&self.additional_data)?;
        Ok(())
    }

    /// Reads params from a buffer.
    pub fn read<R: io::Read>(mut reader: R) -> io::Result<Self> {
        let mut k = [0u8; 4];
        reader.read_exact(&mut k[..])?;
        let k = u32::from_le_bytes(k);
        let n = 1 << k;

        let g: Vec<C> = (0..n)
            .map(|_| C::read(&mut reader))
            .collect::<Result<_, _>>()?;
        let g_lagrange: Vec<C> = (0..n)
            .map(|_| C::read(&mut reader))
            .collect::<Result<_, _>>()?;

        let mut additional_data_len = [0u8; 4];
        reader.read_exact(&mut additional_data_len[..])?;
        let additional_data_len = u32::from_le_bytes(additional_data_len);
        let mut additional_data = vec![0u8; additional_data_len as usize];

        reader.read_exact(&mut additional_data[..])?;

        Ok(Params {
            k,
            n,
            g,
            g_lagrange,
            additional_data,
        })
    }

    /// Returns verifier params with size of Lagrange bases equal to number of public inputs
    pub fn verifier<E: Engine<G1Affine = C>>(
        &self,
        public_inputs_size: usize,
    ) -> io::Result<ParamsVerifier<E>> {
        assert!(public_inputs_size < self.n as usize);
        let g_lagrange = self.g_lagrange[..public_inputs_size].to_vec();
        let g2 = <E::G2Affine as PrimeCurveAffine>::generator();

        let additional_data = self.additional_data.clone();

        let s_g2 = E::G2Affine::read(&mut additional_data.as_slice())?;

        Ok(ParamsVerifier {
            k: self.k,
            n: self.n,
            g1: self.g[0],
            g_lagrange,
            g2,
            s_g2,
        })
    }
}

/// Wrapper type around a blinding factor.
#[derive(Copy, Clone, Eq, PartialEq, Debug)]
pub struct Blind<F>(pub F);

impl<F: FieldExt> Default for Blind<F> {
    fn default() -> Self {
        Blind(F::one())
    }
}

impl<F: FieldExt> Add for Blind<F> {
    type Output = Self;

    fn add(self, rhs: Blind<F>) -> Self {
        Blind(self.0 + rhs.0)
    }
}

impl<F: FieldExt> Mul for Blind<F> {
    type Output = Self;

    fn mul(self, rhs: Blind<F>) -> Self {
        Blind(self.0 * rhs.0)
    }
}

impl<F: FieldExt> AddAssign for Blind<F> {
    fn add_assign(&mut self, rhs: Blind<F>) {
        self.0 += rhs.0;
    }
}

impl<F: FieldExt> MulAssign for Blind<F> {
    fn mul_assign(&mut self, rhs: Blind<F>) {
        self.0 *= rhs.0;
    }
}

impl<F: FieldExt> AddAssign<F> for Blind<F> {
    fn add_assign(&mut self, rhs: F) {
        self.0 += rhs;
    }
}

impl<F: FieldExt> MulAssign<F> for Blind<F> {
    fn mul_assign(&mut self, rhs: F) {
        self.0 *= rhs;
    }
}

impl<E: Engine> ParamsVerifier<E> {
    /// Returns maximum public input size allowed
    pub fn public_inputs_size(&self) -> usize {
        self.g_lagrange.len()
    }

    /// Generates an empty multiscalar multiplication struct using the
    /// appropriate params.
    pub fn empty_msm(&self) -> MSM<E::G1Affine> {
        MSM::new()
    }

    /// Commits to a polynomial using its evaluations over the $2^k$ size
    /// evaluation domain.
    pub fn commit_lagrange(&self, scalars: Vec<E::Scalar>) -> E::G1 {
        let bases = &self.g_lagrange;
        let size = scalars.len();
        assert!(bases.len() >= size);
        best_multiexp(&scalars, &bases[0..size])
    }

    /// Writes params to a buffer.
    pub fn write<W: io::Write>(&self, writer: &mut W) -> io::Result<()> {
        writer.write_all(&self.k.to_le_bytes())?;
        let public_inputs_size = self.public_inputs_size() as u32;
        writer.write_all(&public_inputs_size.to_le_bytes())?;

        writer.write_all(self.g1.to_bytes().as_ref())?;
        writer.write_all(self.g2.to_bytes().as_ref())?;
        writer.write_all(self.s_g2.to_bytes().as_ref())?;
        for el in &self.g_lagrange {
            writer.write_all(el.to_bytes().as_ref())?;
        }
        Ok(())
    }

    /// Reads params from a buffer.
    pub fn read<R: io::Read>(mut reader: R) -> io::Result<Self> {
        let mut k = [0u8; 4];
        reader.read_exact(&mut k[..])?;
        let k = u32::from_le_bytes(k);

        let mut public_inputs_size = [0u8; 4];
        reader.read_exact(&mut public_inputs_size[..])?;
        let public_inputs_size = u32::from_le_bytes(public_inputs_size);

        let n = 1 << k;

        let g1 = E::G1Affine::read(&mut reader)?;
        let g2 = E::G2Affine::read(&mut reader)?;
        let s_g2 = E::G2Affine::read(&mut reader)?;
        let g_lagrange: Vec<E::G1Affine> = (0..public_inputs_size)
            .map(|_| E::G1Affine::read(&mut reader))
            .collect::<Result<_, _>>()?;

        Ok(ParamsVerifier {
            k,
            n,
            g1,
            g2,
            s_g2,
            g_lagrange,
        })
    }
}

/// This will use multithreading if beneficial.
pub fn prev_fft<G: Group>(a: &mut [G], omega: G::Scalar, log_n: u32) {
    let threads = current_num_threads();
    let log_threads = log2_floor(threads);

    if log_n <= log_threads {
        serial_fft(a, omega, log_n);
    } else {
        parallel_fft(a, omega, log_n, log_threads);
    }
}

fn serial_fft<G: Group>(a: &mut [G], omega: G::Scalar, log_n: u32) {
    fn bitreverse(mut n: u32, l: u32) -> u32 {
        let mut r = 0;
        for _ in 0..l {
            r = (r << 1) | (n & 1);
            n >>= 1;
        }
        r
    }

    let n = a.len() as u32;
    assert_eq!(n, 1 << log_n);

    for k in 0..n {
        let rk = bitreverse(k, log_n);
        if k < rk {
            a.swap(rk as usize, k as usize);
        }
    }

    let mut m = 1;
    for _ in 0..log_n {
        let w_m = omega.pow_vartime(&[u64::from(n / (2 * m)), 0, 0, 0]);

        let mut k = 0;
        while k < n {
            let mut w = G::Scalar::one();
            for j in 0..m {
                let mut t = a[(k + j + m) as usize];
                t.group_scale(&w);
                a[(k + j + m) as usize] = a[(k + j) as usize];
                a[(k + j + m) as usize].group_sub(&t);
                a[(k + j) as usize].group_add(&t);
                w *= &w_m;
            }

            k += 2 * m;
        }

        m *= 2;
    }
}

fn parallel_fft<G: Group>(a: &mut [G], omega: G::Scalar, log_n: u32, log_threads: u32) {
    assert!(log_n >= log_threads);

    let num_threads = 1 << log_threads;
    let log_new_n = log_n - log_threads;
    let mut tmp = vec![vec![G::group_zero(); 1 << log_new_n]; num_threads];
    let new_omega = omega.pow_vartime(&[num_threads as u64, 0, 0, 0]);

    scope(|scope| {
        let a = &*a;

        for (j, tmp) in tmp.iter_mut().enumerate() {
            scope.spawn(move |_| {
                // Shuffle into a sub-FFT
                let omega_j = omega.pow_vartime(&[j as u64, 0, 0, 0]);
                let omega_step = omega.pow_vartime(&[(j as u64) << log_new_n, 0, 0, 0]);

                let mut elt = G::Scalar::one();

                for (i, tmp) in tmp.iter_mut().enumerate() {
                    for s in 0..num_threads {
                        let idx = (i + (s << log_new_n)) % (1 << log_n);
                        let mut t = a[idx];
                        t.group_scale(&elt);
                        tmp.group_add(&t);
                        elt *= &omega_step;
                    }
                    elt *= &omega_j;
                }

                // Perform sub-FFT
                serial_fft(tmp, new_omega, log_new_n);
            });
        }
    });

    // Unshuffle
    let mask = (1 << log_threads) - 1;
    for (idx, a) in a.iter_mut().enumerate() {
        *a = tmp[idx & mask][idx >> log_threads];
    }
}

fn log2_floor(num: usize) -> u32 {
    assert!(num > 0);

    let mut pow = 0;

    while (1 << (pow + 1)) <= num {
        pow += 1;
    }

    pow
}

#[cfg(test)]
use pairing::bn256::{Bn256, Fr, G1Affine};

#[test]
fn test_parameter_serialization() {
    const K: u32 = 4;

    let params0: Params<G1Affine> = Params::<G1Affine>::unsafe_setup::<Bn256>(K);

    let mut data: Vec<u8> = Vec::new();
    params0.write(&mut data).unwrap();
    let params1: Params<G1Affine> = Params::read(&data[..]).unwrap();

    assert_eq!(params0.k, params1.k);
    assert_eq!(params0.n, params1.n);
    assert_eq!(params0.g.len(), params1.g.len());
    assert_eq!(params0.g_lagrange.len(), params1.g_lagrange.len());

    assert_eq!(params0.g, params1.g);
    assert_eq!(params0.g_lagrange, params1.g_lagrange);
    assert_eq!(params0.additional_data, params1.additional_data);

    let public_inputs_size = 2;
    let verifier_params0: ParamsVerifier<Bn256> = params0.verifier(public_inputs_size).unwrap();

    assert_eq!(verifier_params0.k, params1.k);
    assert_eq!(verifier_params0.n, params1.n);
    assert_eq!(verifier_params0.g_lagrange.len(), public_inputs_size);
    assert_eq!(
        verifier_params0.s_g2.to_bytes().as_ref(),
        params0.additional_data
    );

    let mut data: Vec<u8> = Vec::new();
    verifier_params0.write(&mut data).unwrap();
    let verifier_params1: ParamsVerifier<Bn256> = ParamsVerifier::read(&data[..]).unwrap();
    assert_eq!(verifier_params0.k, verifier_params1.k);
    assert_eq!(verifier_params0.n, verifier_params1.n);
    assert_eq!(verifier_params0.g1, verifier_params1.g1);
    assert_eq!(verifier_params0.g2, verifier_params1.g2);
    assert_eq!(verifier_params0.s_g2, verifier_params1.s_g2);
    assert_eq!(verifier_params0.g_lagrange, verifier_params1.g_lagrange);
}

#[test]
fn test_commit_lagrange() {
    const K: u32 = 6;

    let params: Params<G1Affine> = Params::<G1Affine>::unsafe_setup::<Bn256>(K);
    let domain = super::EvaluationDomain::new(1, K);

    let mut a = domain.empty_lagrange();

    for (i, a) in a.iter_mut().enumerate() {
        *a = Fr::from(i as u64);
    }

    let b = domain.lagrange_to_coeff(a.clone());
    assert_eq!(params.commit(&b), params.commit_lagrange(&a));
}
