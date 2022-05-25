//! This module contains an implementation of the polynomial commitment scheme
//! described in the [Halo][halo] paper.
//!
//! [halo]: https://eprint.iacr.org/2019/1021

use crate::arithmetic::{
    best_fft, best_multiexp, parallelize, CurveAffine, CurveExt, FieldExt, Group,
};
use crate::helpers::CurveRead;
use crate::poly::commitment::{Blind, CommitmentScheme, Params, ParamsProver, ParamsVerifier, MSM};
use crate::poly::{Coeff, LagrangeCoeff, Polynomial};

use ff::{Field, PrimeField};
use group::{prime::PrimeCurveAffine, Curve, Group as _};
use halo2curves::pairing::Engine;
use std::fmt::Debug;
use std::marker::PhantomData;
use std::ops::{Add, AddAssign, Mul, MulAssign};

use std::io;

use super::msm::MSMKZG;

/// These are the public parameters for the polynomial commitment scheme.
#[derive(Debug, Clone)]
pub struct ParamsKZG<E: Engine> {
    pub(crate) k: u32,
    pub(crate) n: u64,
    pub(crate) g: Vec<E::G1Affine>,
    pub(crate) g_lagrange: Vec<E::G1Affine>,
    pub(crate) g2: E::G2Affine,
    pub(crate) s_g2: E::G2Affine,
}

/// Umbrella commitment scheme construction for all KZG variants
#[derive(Debug)]
pub struct KZGCommitmentScheme<E: Engine> {
    _marker: PhantomData<E>,
}

impl<'params, E: Engine + Debug> CommitmentScheme<'params> for KZGCommitmentScheme<E> {
    type Scalar = E::Scalar;
    type Curve = E::G1Affine;
    type MSM = MSMKZG<E>;

    type ParamsProver = ParamsKZG<E>;
    type ParamsVerifier = ParamsVerifierKZG<E>;

    fn new_params(k: u32) -> Self::ParamsProver {
        ParamsKZG::new(k)
    }

    fn read_params<R: io::Read>(reader: &mut R) -> io::Result<Self::ParamsProver> {
        ParamsKZG::read(reader)
    }
}

// TODO: see the issue at https://github.com/appliedzkp/halo2/issues/45
// So we probably need much smaller verifier key. However for new bases in g1 should be in verifier keys.
/// KZG multi-open verification parameters
pub type ParamsVerifierKZG<C> = ParamsKZG<C>;

impl<'params, E: Engine + Debug> Params<'params, E::G1Affine, MSMKZG<E>> for ParamsKZG<E> {
    fn k(&self) -> u32 {
        self.k
    }

    fn n(&self) -> u64 {
        self.n
    }

    fn empty_msm(&'params self) -> MSMKZG<E> {
        MSMKZG::new()
    }

    fn commit_lagrange(
        &self,
        poly: &Polynomial<E::Scalar, LagrangeCoeff>,
        _: Blind<E::Scalar>,
    ) -> E::G1 {
        let mut scalars = Vec::with_capacity(poly.len());
        scalars.extend(poly.iter());
        let bases = &self.g_lagrange;
        let size = scalars.len();
        assert!(bases.len() >= size);
        best_multiexp(&scalars, &bases[0..size])
    }

    /// Writes params to a buffer.
    fn write<W: io::Write>(&self, writer: &mut W) -> io::Result<()> {
        use group::GroupEncoding;
        writer.write_all(&self.k.to_le_bytes())?;
        for el in self.g.iter() {
            writer.write_all(el.to_bytes().as_ref())?;
        }
        for el in self.g_lagrange.iter() {
            writer.write_all(el.to_bytes().as_ref())?;
        }
        writer.write_all(self.g2.to_bytes().as_ref())?;
        writer.write_all(self.s_g2.to_bytes().as_ref())?;
        Ok(())
    }

    /// Reads params from a buffer.
    fn read<R: io::Read>(reader: &mut R) -> io::Result<Self> {
        use group::GroupEncoding;

        let mut k = [0u8; 4];
        reader.read_exact(&mut k[..])?;
        let k = u32::from_le_bytes(k);
        let n = 1 << k;

        let load_points_from_file_parallelly = |reader: &mut R| -> io::Result<Vec<E::G1Affine>> {
            let mut points_compressed =
                vec![<<E as Engine>::G1Affine as GroupEncoding>::Repr::default(); n];
            for points_compressed in points_compressed.iter_mut() {
                reader.read_exact((*points_compressed).as_mut())?;
            }

            let points = points_compressed
                .iter()
                .map(|compressed| Option::from(E::G1Affine::from_bytes(compressed)).unwrap())
                .collect();

            Ok(points)
        };

        let g = load_points_from_file_parallelly(reader)?;
        let g_lagrange = load_points_from_file_parallelly(reader)?;

        let mut point_compressed = <<E as Engine>::G2Affine as GroupEncoding>::Repr::default();
        reader.read_exact(point_compressed.as_mut())?;
        let g2 = Option::from(E::G2Affine::from_bytes(&point_compressed)).unwrap();

        let mut point_compressed = <<E as Engine>::G2Affine as GroupEncoding>::Repr::default();
        reader.read_exact(point_compressed.as_mut())?;
        let s_g2 = Option::from(E::G2Affine::from_bytes(&point_compressed)).unwrap();

        Ok(Self {
            k,
            n: n as u64,
            g,
            g_lagrange,
            g2,
            s_g2,
        })
    }
}

impl<'params, E: Engine + Debug> ParamsVerifier<'params, E::G1Affine, MSMKZG<E>> for ParamsKZG<E> {}

impl<'params, E: Engine + Debug> ParamsProver<'params, E::G1Affine, MSMKZG<E>> for ParamsKZG<E> {
    // TODO: Verifier requires shorter key. Keeping as same as prover params since PLONK
    // implementation here requires commitment to public inputs.
    type ParamsVerifier = ParamsVerifierKZG<E>;

    fn verifier_params(&'params self) -> &'params Self::ParamsVerifier {
        self
    }

    fn new(k: u32) -> Self {
        // Largest root of unity exponent of the Engine is `2^E::Scalar::S`, so we can
        // only support FFTs of polynomials below degree `2^E::Scalar::S`.
        assert!(k <= E::Scalar::S);
        let n: u64 = 1 << k;

        // Calculate g = [G1, [s] G1, [s^2] G1, ..., [s^(n-1)] G1] in parallel.
        let g1 = E::G1Affine::generator();
        use rand_core::OsRng;
        let s = <E::Scalar>::random(OsRng);

        let mut g_projective = vec![E::G1::group_zero(); n as usize];
        parallelize(&mut g_projective, |g, start| {
            let mut current_g: E::G1 = g1.into();
            current_g *= s.pow_vartime(&[start as u64]);
            for g in g.iter_mut() {
                *g = current_g;
                current_g *= s;
            }
        });

        let g = {
            let mut g = vec![E::G1Affine::identity(); n as usize];
            parallelize(&mut g, |g, starts| {
                E::G1::batch_normalize(&g_projective[starts..(starts + g.len())], g);
            });
            g
        };

        let mut g_lagrange_projective = vec![E::G1::group_zero(); n as usize];
        let mut root = E::Scalar::ROOT_OF_UNITY_INV.invert().unwrap();
        for _ in k..E::Scalar::S {
            root = root.square();
        }
        let n_inv = Option::<E::Scalar>::from(E::Scalar::from(n).invert())
            .expect("inversion should be ok for n = 1<<k");
        let multiplier = (s.pow_vartime(&[n as u64]) - E::Scalar::one()) * n_inv;
        parallelize(&mut g_lagrange_projective, |g, start| {
            for (idx, g) in g.iter_mut().enumerate() {
                let offset = start + idx;
                let root_pow = root.pow_vartime(&[offset as u64]);
                let scalar = multiplier * root_pow * (s - root_pow).invert().unwrap();
                *g = g1 * scalar;
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
        let s_g2 = (g2 * s).into();

        ParamsKZG {
            k,
            n,
            g,
            g_lagrange,
            g2,
            s_g2,
        }
    }

    fn commit(&self, poly: &Polynomial<E::Scalar, Coeff>, _: Blind<E::Scalar>) -> E::G1 {
        let mut scalars = Vec::with_capacity(poly.len());
        scalars.extend(poly.iter());
        let bases = &self.g;
        let size = scalars.len();
        assert!(bases.len() >= size);
        best_multiexp(&scalars, &bases[0..size])
    }

    fn get_g(&self) -> Vec<E::G1Affine> {
        self.g.clone()
    }
}

#[cfg(test)]
mod test {
    // TODO: duplicate is in IPA side

    use crate::arithmetic::{
        best_fft, best_multiexp, parallelize, CurveAffine, CurveExt, FieldExt, Group,
    };
    use crate::helpers::CurveRead;
    use crate::poly::commitment::ParamsProver;
    use crate::poly::commitment::{Blind, CommitmentScheme, Params, MSM};
    use crate::poly::kzg::commitment::ParamsKZG;
    use crate::poly::kzg::msm::MSMKZG;
    use crate::poly::kzg::multiopen::ProverSHPLONK;
    use crate::poly::{Coeff, LagrangeCoeff, Polynomial};

    use ff::{Field, PrimeField};
    use group::{prime::PrimeCurveAffine, Curve, Group as _};
    use halo2curves::bn256::G1Affine;
    use std::marker::PhantomData;
    use std::ops::{Add, AddAssign, Mul, MulAssign};

    use std::io;

    #[test]
    fn test_commit_lagrange() {
        const K: u32 = 6;

        use rand_core::OsRng;

        use crate::poly::EvaluationDomain;
        use halo2curves::bn256::{Bn256, Fr};

        let params = ParamsKZG::<Bn256>::new(K);
        let domain = EvaluationDomain::new(1, K);

        let mut a = domain.empty_lagrange();

        for (i, a) in a.iter_mut().enumerate() {
            *a = Fr::from(i as u64);
        }

        let b = domain.lagrange_to_coeff(a.clone());

        let alpha = Blind(Fr::random(OsRng));

        assert_eq!(params.commit(&b, alpha), params.commit_lagrange(&a, alpha));
    }
}
