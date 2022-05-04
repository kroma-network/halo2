use crate::arithmetic::{parallelize, CurveAffine, Engine};
use group::prime::PrimeCurveAffine;
use rand_core::OsRng;

/// These are the prover parameters for the polynomial commitment scheme.
#[derive(Debug)]
pub struct ProverParams<C: CurveAffine> {
    // polynomial degree log n
    pub(crate) k: u32,
    // polynomial degree n
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

impl<C: CurveAffine> ProverParams<C> {
    pub fn unsafe_setup<E: Engine>(k: u32) -> ProverParams<E::G1Affine> {
        // Largest root of unity exponent of the Engine is `2^E::Scalar::S`, so we can
        // only support FFTs of polynomials below degree `2^E::Scalar::S`.
        assert!(k <= E::Scalar::S);
        let n: u64 = 1 << k;
        let g1 = <E::G1Affine as PrimeCurveAffine>::generator();

        // Generate randomness `s` for srs
        let s = E::Scalar::random(OsRng);

        // Calculate g = [G1, [s] G1, [s^2] G1, ..., [s^(n-1)] G1] in parallel.
        // projective form of coordinate
        let mut g_projective = vec![E::G1::group_zero(); n as usize];
        parallelize(&mut g_projective, |g, start| {
            let mut current_g: E::G1 = g1.into();
            current_g *= s.pow_vartime(&[start as u64]);
            for g in g.iter_mut() {
                *g = current_g;
                current_g *= s;
            }
        });

        // Convert g_projective elements into affine elements `g`
        // affine form of coordinate
        let mut g = vec![E::G1Affine::identity(); n as usize];
        parallelize(&mut g, |g, starts| {
            E::G1::batch_normalize(&g_projective[starts..(starts + g.len())], g);
        });

        // Calculate Lagrange form srs
        let mut g_lagrange_projective = vec![E::G1::group_zero(); n as usize];
        let mut root = E::Scalar::ROOT_OF_UNITY_INV.invert().unwrap();
        for _ in k..E::Scalar::S {
            root = root.square();
        }
        let n_inv = Option::<E::Scalar>::from(E::Scalar::from(n).invert())
            .expect("inversion should be ok for n = 1<<k");
        let multiplier = (s.pow_vartime(&[n as u64]) - E::Scalar::one()) * n_inv;

        // Evaluate all of the Lagrange basis polynomials using an inverse FFT.
        parallelize(&mut g_lagrange_projective, |g, start| {
            for (idx, g) in g.iter_mut().enumerate() {
                let offset = start + idx;
                let root_pow = root.pow_vartime(&[offset as u64]);
                let scalar = multiplier * root_pow * (s - root_pow).invert().unwrap();
                *g = g1 * scalar;
            }
        });

        // Convert g_lagrange_projective elements into affine elements `g_lagrange`
        // affine form of coordinate
        let mut g_lagrange = vec![E::G1Affine::identity(); n as usize];
        parallelize(&mut g_lagrange, |g_lagrange, starts| {
            E::G1::batch_normalize(
                &g_lagrange_projective[starts..(starts + g_lagrange.len())],
                g_lagrange,
            );
        });

        let g2 = <E::G2Affine as PrimeCurveAffine>::generator();
        let s_g2 = g2 * s;
        let additional_data = Vec::from(s_g2.to_bytes().as_ref());
        Self {
            k,
            n,
            g,
            g_lagrange,
            additional_data,
        }
    }

    /// This computes a commitment to a polynomial described by the provided
    /// slice of coefficients.
    pub fn commit(&self, poly: &Polynomial<C::Scalar, Coeff>) -> C::Curve {
        let mut scalars = Vec::with_capacity(poly.len());
        scalars.extend(poly.iter());
        let bases = &self.g;
        let size = scalars.len();
        assert!(bases.len() >= size);
        best_multiexp(&scalars, &bases[0..size])
    }

    /// This commits to a polynomial using its evaluations over the $2^k$ size
    /// evaluation domain.
    pub fn commit_lagrange(&self, poly: &Polynomial<C::Scalar, LagrangeCoeff>) -> C::Curve {
        let mut scalars = Vec::with_capacity(poly.len());
        scalars.extend(poly.iter());
        let bases = &self.g_lagrange;
        let size = scalars.len();
        assert!(bases.len() >= size);
        best_multiexp(&scalars, &bases[0..size])
    }
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
