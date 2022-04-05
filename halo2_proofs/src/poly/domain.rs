//! Contains utilities for performing polynomial arithmetic over an evaluation
//! domain that is of a suitable size for the application.

use crate::{
    arithmetic::{best_fft, parallelize, FieldExt, Group},
    multicore,
    plonk::Assigned,
};

use super::{Coeff, ExtendedLagrangeCoeff, LagrangeCoeff, Polynomial, Rotation};

use group::ff::{BatchInvert, Field, PrimeField};

use std::{marker::PhantomData, os::unix::prelude::FileExt};

/// This structure hold the twiddles and radix for each layer
#[derive(Debug)]
pub struct FFTBitReverseCache {
    /// indexes for bit reverse
    pub a_indexes: Vec<usize>,
    /// indexes for bit reverse
    pub b_indexes: Vec<usize>,
}

impl FFTBitReverseCache {
    /// bit reverse method
    pub fn sort_bit_reverse<F: FieldExt>(&self, a: &mut [F]) {
        multicore::scope(|s| {
            s.spawn(move |_| {
                for (f, s) in self.a_indexes.iter().zip(self.b_indexes.iter()) {
                    a.swap(*f, *s);
                }
            })
        })
    }
}

/// This structure hold the twiddles and radix for each layer
#[derive(Debug)]
pub struct FFTButterlyCache<F: FieldExt> {
    /// n half
    pub half: usize,
    /// stages
    pub layer: usize,
    /// twiddles
    pub twiddles: Vec<F>,
    /// odd k flag
    pub is_odd: bool,
    /// low degree flag
    pub is_low: bool,
}

impl<F: FieldExt> FFTButterlyCache<F> {
    /// fft
    pub fn butterfly_arithmetic(&self, a: &mut [F], n: usize) {
        if self.is_low {
            for i in 0..n / 2 {
                let t = a[2 * i + 1];
                a[2 * i + 1] = a[2 * i];
                a[2 * i] += t;
                a[2 * i + 1] -= t;
            }

            match n {
                4 => {
                    let tw = a[2];
                    a[2] = a[0];
                    a[0] += tw;
                    a[2] -= tw;
                    let tw = a[3] * self.twiddles[1];
                    a[3] = a[1];
                    a[1] += tw;
                    a[3] -= tw;
                }
                8 => {
                    let tw = a[2];
                    a[2] = a[0];
                    a[0] += tw;
                    a[2] -= tw;
                    let tw = a[3] * self.twiddles[2];
                    a[3] = a[1];
                    a[1] += tw;
                    a[3] -= tw;
                    let tw = a[6];
                    a[6] = a[4];
                    a[4] += tw;
                    a[6] -= tw;
                    let tw = a[7] * self.twiddles[2];
                    a[7] = a[5];
                    a[5] += tw;
                    a[7] -= tw;

                    let tw = a[4];
                    a[4] = a[0];
                    a[0] += tw;
                    a[4] -= tw;
                    let tw = a[5] * self.twiddles[1];
                    a[5] = a[1];
                    a[1] += tw;
                    a[5] -= tw;
                    let tw = a[6] * self.twiddles[2];
                    a[6] = a[2];
                    a[2] += tw;
                    a[6] -= tw;
                    let tw = a[7] * self.twiddles[3];
                    a[7] = a[3];
                    a[3] += tw;
                    a[7] -= tw;
                }
                _ => {}
            }
        } else {
            let mut elements = if self.is_odd { 16 } else { 32 };

            // twiddles factor 16th root of unity
            let tw_offset = self.half / 8;
            let tw_1 = self.twiddles[tw_offset];
            let tw_2 = self.twiddles[tw_offset * 2];
            let tw_3 = self.twiddles[tw_offset * 3];
            let tw_4 = self.twiddles[tw_offset * 4];
            let tw_5 = self.twiddles[tw_offset * 5];
            let tw_6 = self.twiddles[tw_offset * 6];
            let tw_7 = self.twiddles[tw_offset * 7];

            for a in a.chunks_mut(16) {
                // first layer butterfly arithmetic
                let a_a = a[0] + a[1];
                let a_b = a[0] - a[1];
                let a_c = a[2] + a[3];
                let a_d = a[2] - a[3];
                let b_a = a[4] + a[5];
                let b_b = a[4] - a[5];
                let b_c = a[6] + a[7];
                let b_d = a[6] - a[7];
                let c_a = a[8] + a[9];
                let c_b = a[8] - a[9];
                let c_c = a[10] + a[11];
                let c_d = a[10] - a[11];
                let d_a = a[12] + a[13];
                let d_b = a[12] - a[13];
                let d_c = a[14] + a[15];
                let d_d = a[14] - a[15];

                // second layer butterfly airthmetic
                let a_w_0 = a_d * tw_4;
                let b_w_0 = b_d * tw_4;
                let c_w_0 = c_d * tw_4;
                let d_w_0 = d_d * tw_4;

                let e_a = a_a + a_c;
                let e_b = a_b + a_w_0;
                let e_c = a_a - a_c;
                let e_d = a_b - a_w_0;
                let f_a = b_a + b_c;
                let f_b = b_b + b_w_0;
                let f_c = b_a - b_c;
                let f_d = b_b - b_w_0;
                let g_a = c_a + c_c;
                let g_b = c_b + c_w_0;
                let g_c = c_a - c_c;
                let g_d = c_b - c_w_0;
                let h_a = d_a + d_c;
                let h_b = d_b + d_w_0;
                let h_c = d_a - d_c;
                let h_d = d_b - d_w_0;

                // third layer butterfly airthmetic
                let f_w_1 = f_b * tw_2;
                let f_w_2 = f_c * tw_4;
                let f_w_3 = f_d * tw_6;
                let h_w_1 = h_b * tw_2;
                let h_w_2 = h_c * tw_4;
                let h_w_3 = h_d * tw_6;

                let i_a = e_a + f_a;
                let i_b = e_b + f_w_1;
                let i_c = e_c + f_w_2;
                let i_d = e_d + f_w_3;
                let j_a = e_a - f_a;
                let j_b = e_b - f_w_1;
                let j_c = e_c - f_w_2;
                let j_d = e_d - f_w_3;
                let k_a = g_a + h_a;
                let k_b = g_b + h_w_1;
                let k_c = g_c + h_w_2;
                let k_d = g_d + h_w_3;
                let l_a = g_a - h_a;
                let l_b = g_b - h_w_1;
                let l_c = g_c - h_w_2;
                let l_d = g_d - h_w_3;

                if self.is_odd {
                    a[0] = i_a;
                    a[1] = i_b;
                    a[2] = i_c;
                    a[3] = i_d;
                    a[4] = j_a;
                    a[5] = j_b;
                    a[6] = j_c;
                    a[7] = j_d;
                    a[8] = k_a;
                    a[9] = k_b;
                    a[10] = k_c;
                    a[11] = k_d;
                    a[12] = l_a;
                    a[13] = l_b;
                    a[14] = l_c;
                    a[15] = l_d;
                } else {
                    // forth layer butterfly airthmetic
                    let k_w_1 = k_b * tw_1;
                    let k_w_2 = k_c * tw_2;
                    let k_w_3 = k_d * tw_3;
                    let k_w_4 = l_a * tw_4;
                    let k_w_5 = l_b * tw_5;
                    let k_w_6 = l_c * tw_6;
                    let k_w_7 = l_d * tw_7;

                    a[0] = i_a + k_a;
                    a[1] = i_b + k_w_1;
                    a[2] = i_c + k_w_2;
                    a[3] = i_d + k_w_3;
                    a[4] = j_a + k_w_4;
                    a[5] = j_b + k_w_5;
                    a[6] = j_c + k_w_6;
                    a[7] = j_d + k_w_7;
                    a[8] = i_a - k_a;
                    a[9] = i_b - k_w_1;
                    a[10] = i_c - k_w_2;
                    a[11] = i_d - k_w_3;
                    a[12] = j_a - k_w_4;
                    a[13] = j_b - k_w_5;
                    a[14] = j_c - k_w_6;
                    a[15] = j_d - k_w_7;
                }
            }

            // four radix butterfly arithmetic
            for _ in 0..self.layer {
                // element and twiddles offset
                let offset = elements / 2;
                let tw_offset = self.half / elements;
                multicore::scope(|scope| {
                    for input in a.chunks_mut(elements * 2) {
                        scope.spawn(move |_| {
                            for p in 0..offset {
                                // indexes of this layer
                                let first = p;
                                let second = first + offset;
                                let third = second + offset;
                                let fourth = third + offset;

                                // twiddle factor arithmetic for upper side
                                let a_tw_idx = tw_offset * 2 * p;
                                let second_tw = input[second] * self.twiddles[a_tw_idx];
                                let fourth_tw = input[fourth] * self.twiddles[a_tw_idx];

                                // upper side butterfly arithmetic
                                let a = input[first] + second_tw;
                                let b = input[first] - second_tw;
                                let c = input[third] + fourth_tw;
                                let d = input[third] - fourth_tw;

                                // twiddle factor arithmetic for bottom side
                                let b_tw1_idx = a_tw_idx / 2;
                                let b_tw2_idx = b_tw1_idx + self.half / 2;
                                let c_tw = c * self.twiddles[b_tw1_idx];
                                let d_tw = d * self.twiddles[b_tw2_idx];

                                // bottom side butterfly arithmetic
                                input[first] = a + c_tw;
                                input[second] = b + d_tw;
                                input[third] = a - c_tw;
                                input[fourth] = b - d_tw;
                            }
                        });
                    }
                });
                elements *= 4;
            }
        }
    }
}

/// This structure hold the twiddles and radix for each layer
#[derive(Debug)]
pub struct FFTCache<F: FieldExt> {
    /// n half
    pub bit_reverse: FFTBitReverseCache,
    /// stages
    pub butterfly: FFTButterlyCache<F>,
}

/// This structure hold the twiddles and radix for each layer
#[derive(Debug)]
pub struct FFTHelper<F: FieldExt> {
    /// fft data
    pub fft_data: FFTCache<F>,
    /// extended fft data
    pub ext_fft_data: FFTCache<F>,
}

impl<F: FieldExt> FFTHelper<F> {
    /// `FFTHelper` init
    pub fn new(k: usize, ext_k: usize, omega: F, ext_omega: F) -> Self {
        fn bitreverse(mut n: usize, l: usize) -> usize {
            let mut r = 0;
            for _ in 0..l {
                r = (r << 1) | (n & 1);
                n >>= 1;
            }
            r
        }

        fn get_fft_helper(k: usize) -> (usize, usize, bool, usize, bool) {
            let n = 1 << k;
            let half = n / 2;
            let (offset, is_low) = if k < 4 { (0, true) } else { (k - 4, false) };
            let is_odd = k % 2 == 1;
            let layer = if offset == 0 {
                0
            } else {
                offset / 2 + if is_odd { 1 } else { 0 }
            };
            (n, half, is_low, layer, is_odd)
        }

        let (n, half, is_low, layer, is_odd) = get_fft_helper(k);
        let (ext_n, ext_half, ext_is_low, ext_layer, ext_is_odd) = get_fft_helper(ext_k);

        // calculate twiddles factor
        let mut w = F::one();
        let mut e_w = F::one();
        let mut twiddles = vec![F::one(); half];
        let mut ext_twiddles = vec![F::one(); ext_half];

        // init bit reverse indexes
        let mut a_indexes = Vec::with_capacity(half / 2);
        let mut b_indexes = Vec::with_capacity(half / 2);
        let mut ext_a_indexes = Vec::with_capacity(ext_half / 2);
        let mut ext_b_indexes = Vec::with_capacity(ext_half / 2);

        // mutable reference for multicore
        let stash = &mut twiddles;
        let e_stash = &mut ext_twiddles;
        let tmp_a_idx = &mut a_indexes;
        let tmp_b_idx = &mut b_indexes;
        let tmp_ext_a_idx = &mut ext_a_indexes;
        let tmp_ext_b_idx = &mut ext_b_indexes;

        multicore::scope(|scope| {
            scope.spawn(move |_| {
                for i in 0..n {
                    let ri = bitreverse(i, k);
                    if i < ri {
                        tmp_a_idx.push(ri);
                        tmp_b_idx.push(i);
                    }
                }
            });
            scope.spawn(move |_| {
                for i in 0..ext_n {
                    let ri = bitreverse(i, k);
                    if i < ri {
                        tmp_ext_a_idx.push(ri);
                        tmp_ext_b_idx.push(i);
                    }
                }
            });
            scope.spawn(move |_| {
                for tw in stash.iter_mut() {
                    *tw = w;
                    w *= omega;
                }
            });
            scope.spawn(move |_| {
                for tw in e_stash.iter_mut() {
                    *tw = e_w;
                    e_w *= ext_omega;
                }
            });
        });

        Self {
            fft_data: FFTCache {
                bit_reverse: FFTBitReverseCache {
                    a_indexes,
                    b_indexes,
                },
                butterfly: FFTButterlyCache {
                    half,
                    layer,
                    twiddles,
                    is_odd,
                    is_low,
                },
            },
            ext_fft_data: FFTCache {
                bit_reverse: FFTBitReverseCache {
                    a_indexes: ext_a_indexes,
                    b_indexes: ext_b_indexes,
                },
                butterfly: FFTButterlyCache {
                    half: ext_half,
                    layer: ext_layer,
                    twiddles: ext_twiddles,
                    is_odd: ext_is_odd,
                    is_low: ext_is_low,
                },
            },
        }
    }
}

/// This structure hold the twiddles and radix for each layer
#[derive(Debug)]
pub struct FFTData<F: FieldExt> {
    /// n half
    pub half: usize,
    /// stages
    pub layer: usize,
    /// indexes for bit reverse
    pub indexes: Vec<usize>,
    /// twiddles
    pub twiddles: Vec<F>,
    /// inv twiddles
    pub inv_twiddles: Vec<F>,
    /// odd k flag
    pub is_odd: bool,
    /// low degree flag
    pub is_low: bool,
}

impl<F: FieldExt> FFTData<F> {
    /// Create twiddles and stages data
    pub fn new(n: usize, omega: F, omega_inv: F, k: usize) -> Self {
        let half = n / 2;
        let (offset, is_low) = if k < 4 { (0, true) } else { (k - 4, false) };
        let mut counter = 2;

        // calculate twiddles factor
        let mut w = F::one();
        let mut i = F::one();
        let mut twiddles = vec![F::one(); half];
        let mut inv_twiddles = vec![F::one(); half];

        // init bit reverse indexes
        let mut indexes = vec![0; n];

        // mutable reference for multicore
        let stash = &mut twiddles;
        let i_stash = &mut inv_twiddles;

        multicore::scope(|scope| {
            scope.spawn(move |_| {
                for tw in stash.iter_mut() {
                    *tw = w;
                    w *= omega;
                }
            });
            scope.spawn(move |_| {
                for tw in i_stash.iter_mut() {
                    *tw = i;
                    i *= omega_inv;
                }
            });
        });

        if k % 2 != 0 {
            indexes[0] = 0;
            indexes[1] = 1;
        } else {
            indexes[0] = 0;
            indexes[1] = 2;
            indexes[2] = 1;
            indexes[3] = 3;
            counter *= 2;
        }

        // calculate bit reverse indexes
        while counter != n {
            for i in 0..counter {
                indexes[i] *= 4;
                indexes[i + counter] = indexes[i] + 2;
                indexes[i + 2 * counter] = indexes[i] + 1;
                indexes[i + 3 * counter] = indexes[i] + 3;
            }
            counter *= 4;
        }

        Self {
            half,
            layer: offset / 2,
            indexes,
            twiddles,
            inv_twiddles,
            is_odd: k % 2 == 1,
            is_low,
        }
    }
}

/// This structure hold the twiddles and radix for each layer
#[derive(Debug)]
pub struct EvaluationHelper<F: FieldExt> {
    /// fft data
    pub fft_data: FFTData<F>,
    /// extended fft data
    pub ext_fft_data: FFTData<F>,
}

impl<F: FieldExt> EvaluationHelper<F> {
    /// Create twiddles and stages data
    pub fn new(
        n: usize,
        extended_n: usize,
        omega: F,
        omega_inv: F,
        extended_omega: F,
        extended_omega_inv: F,
        k: usize,
        extended_k: usize,
    ) -> Self {
        assert_eq!(n, 1 << k);
        assert_eq!(extended_n, 1 << extended_k);

        let fft_data = FFTData::new(n, omega, omega_inv, k);
        let ext_fft_data = FFTData::new(extended_n, extended_omega, extended_omega_inv, extended_k);

        Self {
            fft_data,
            ext_fft_data,
        }
    }
}

/// This structure contains precomputed constants and other details needed for
/// performing operations on an evaluation domain of size $2^k$ and an extended
/// domain of size $2^{k} * j$ with $j \neq 0$.
#[derive(Debug)]
pub struct EvaluationDomain<F: FieldExt> {
    n: u64,
    k: u32,
    extended_k: u32,
    omega: F,
    omega_inv: F,
    extended_omega: F,
    extended_omega_inv: F,
    g_coset: F,
    g_coset_inv: F,
    quotient_poly_degree: u64,
    ifft_divisor: F,
    extended_ifft_divisor: F,
    t_evaluations: Vec<F>,
    barycentric_weight: F,
    evaluation_helper: EvaluationHelper<F>,
}

impl<F: FieldExt> EvaluationDomain<F> {
    /// This constructs a new evaluation domain object based on the provided
    /// values $j, k$.
    pub fn new(j: u32, k: u32) -> Self {
        // quotient_poly_degree * params.n - 1 is the degree of the quotient polynomial
        let quotient_poly_degree = (j - 1) as u64;

        // n = 2^k
        let n = 1u64 << k;

        // We need to work within an extended domain, not params.k but params.k + i
        // for some integer i such that 2^(params.k + i) is sufficiently large to
        // describe the quotient polynomial.
        let mut extended_k = k;
        while (1 << extended_k) < (n * quotient_poly_degree) {
            extended_k += 1;
        }

        let mut extended_omega = F::Scalar::root_of_unity();

        // Get extended_omega, the 2^{extended_k}'th root of unity
        // The loop computes extended_omega = omega^{2 ^ (S - extended_k)}
        // Notice that extended_omega ^ {2 ^ extended_k} = omega ^ {2^S} = 1.
        for _ in extended_k..F::Scalar::S {
            extended_omega = extended_omega.square();
        }
        let extended_omega = extended_omega;
        let mut extended_omega_inv = extended_omega; // Inversion computed later

        // Get omega, the 2^{k}'th root of unity (i.e. n'th root of unity)
        // The loop computes omega = extended_omega ^ {2 ^ (extended_k - k)}
        //           = (omega^{2 ^ (S - extended_k)})  ^ {2 ^ (extended_k - k)}
        //           = omega ^ {2 ^ (S - k)}.
        // Notice that omega ^ {2^k} = omega ^ {2^S} = 1.
        let mut omega = extended_omega;
        for _ in k..extended_k {
            omega = omega.square();
        }
        let omega = omega;
        let mut omega_inv = omega; // Inversion computed later

        // We use zeta here because we know it generates a coset, and it's available
        // already.
        // The coset evaluation domain is:
        // zeta {1, extended_omega, extended_omega^2, ..., extended_omega^{(2^extended_k) - 1}}
        let g_coset = F::Scalar::ZETA;
        let g_coset_inv = g_coset.square();

        let mut t_evaluations = Vec::with_capacity(1 << (extended_k - k));
        {
            // Compute the evaluations of t(X) = X^n - 1 in the coset evaluation domain.
            // We don't have to compute all of them, because it will repeat.
            let orig = F::Scalar::ZETA.pow_vartime(&[n as u64, 0, 0, 0]);
            let step = extended_omega.pow_vartime(&[n as u64, 0, 0, 0]);
            let mut cur = orig;
            loop {
                t_evaluations.push(cur);
                cur *= &step;
                if cur == orig {
                    break;
                }
            }
            assert_eq!(t_evaluations.len(), 1 << (extended_k - k));

            // Subtract 1 from each to give us t_evaluations[i] = t(zeta * extended_omega^i)
            for coeff in &mut t_evaluations {
                *coeff -= &F::Scalar::one();
            }

            // Invert, because we're dividing by this polynomial.
            // We invert in a batch, below.
        }

        let mut ifft_divisor = F::Scalar::from(1 << k); // Inversion computed later
        let mut extended_ifft_divisor = F::Scalar::from(1 << extended_k); // Inversion computed later

        // The barycentric weight of 1 over the evaluation domain
        // 1 / \prod_{i != 0} (1 - omega^i)
        let mut barycentric_weight = F::Scalar::from(n); // Inversion computed later

        // Compute batch inversion
        t_evaluations
            .iter_mut()
            .chain(Some(&mut ifft_divisor))
            .chain(Some(&mut extended_ifft_divisor))
            .chain(Some(&mut barycentric_weight))
            .chain(Some(&mut extended_omega_inv))
            .chain(Some(&mut omega_inv))
            .batch_invert();

        EvaluationDomain {
            n,
            k,
            extended_k,
            omega,
            omega_inv,
            extended_omega,
            extended_omega_inv,
            g_coset,
            g_coset_inv,
            quotient_poly_degree,
            ifft_divisor,
            extended_ifft_divisor,
            t_evaluations,
            barycentric_weight,
            evaluation_helper: EvaluationHelper::<F::Scalar>::new(
                n as usize,
                (1 << extended_k) as usize,
                omega,
                omega_inv,
                extended_omega,
                extended_omega_inv,
                k as usize,
                extended_k as usize,
            ),
        }
    }

    /// Obtains a polynomial in Lagrange form when given a vector of Lagrange
    /// coefficients of size `n`; panics if the provided vector is the wrong
    /// length.
    pub fn lagrange_from_vec(&self, values: Vec<F>) -> Polynomial<F, LagrangeCoeff> {
        assert_eq!(values.len(), self.n as usize);

        Polynomial {
            values,
            _marker: PhantomData,
        }
    }

    /// Obtains a polynomial in coefficient form when given a vector of
    /// coefficients of size `n`; panics if the provided vector is the wrong
    /// length.
    pub fn coeff_from_vec(&self, values: Vec<F>) -> Polynomial<F, Coeff> {
        assert_eq!(values.len(), self.n as usize);

        Polynomial {
            values,
            _marker: PhantomData,
        }
    }

    /// Returns an empty (zero) polynomial in the coefficient basis
    pub fn empty_coeff(&self) -> Polynomial<F, Coeff> {
        Polynomial {
            values: vec![F::group_zero(); self.n as usize],
            _marker: PhantomData,
        }
    }

    /// Returns an empty (zero) polynomial in the Lagrange coefficient basis
    pub fn empty_lagrange(&self) -> Polynomial<F, LagrangeCoeff> {
        Polynomial {
            values: vec![F::group_zero(); self.n as usize],
            _marker: PhantomData,
        }
    }

    /// Returns an empty (zero) polynomial in the Lagrange coefficient basis, with
    /// deferred inversions.
    pub(crate) fn empty_lagrange_assigned(&self) -> Polynomial<Assigned<F>, LagrangeCoeff>
    where
        F: Field,
    {
        Polynomial {
            values: vec![F::group_zero().into(); self.n as usize],
            _marker: PhantomData,
        }
    }

    /// Returns a constant polynomial in the Lagrange coefficient basis
    pub fn constant_lagrange(&self, scalar: F) -> Polynomial<F, LagrangeCoeff> {
        Polynomial {
            values: vec![scalar; self.n as usize],
            _marker: PhantomData,
        }
    }

    /// Returns an empty (zero) polynomial in the extended Lagrange coefficient
    /// basis
    pub fn empty_extended(&self) -> Polynomial<F, ExtendedLagrangeCoeff> {
        Polynomial {
            values: vec![F::group_zero(); self.extended_len()],
            _marker: PhantomData,
        }
    }

    /// Returns a constant polynomial in the extended Lagrange coefficient
    /// basis
    pub fn constant_extended(&self, scalar: F) -> Polynomial<F, ExtendedLagrangeCoeff> {
        Polynomial {
            values: vec![scalar; self.extended_len()],
            _marker: PhantomData,
        }
    }

    /// This takes us from an n-length vector into the coefficient form.
    ///
    /// This function will panic if the provided vector is not the correct
    /// length.
    pub fn lagrange_to_coeff(&self, mut a: Polynomial<F, LagrangeCoeff>) -> Polynomial<F, Coeff> {
        assert_eq!(a.values.len(), 1 << self.k);

        // Perform inverse FFT to obtain the polynomial in coefficient form
        Self::ifft(
            &mut a.values,
            &self.evaluation_helper.fft_data,
            true,
            self.ifft_divisor,
        );

        Polynomial {
            values: a.values,
            _marker: PhantomData,
        }
    }

    /// This takes us from an n-length coefficient vector into a coset of the extended
    /// evaluation domain, rotating by `rotation` if desired.
    pub fn coeff_to_extended(
        &self,
        mut a: Polynomial<F, Coeff>,
    ) -> Polynomial<F, ExtendedLagrangeCoeff> {
        assert_eq!(a.values.len(), 1 << self.k);

        self.distribute_powers_zeta(&mut a.values, true);
        a.values.resize(self.extended_len(), F::group_zero());
        best_fft(&mut a.values, &self.evaluation_helper.ext_fft_data, false);

        Polynomial {
            values: a.values,
            _marker: PhantomData,
        }
    }

    /// Rotate the extended domain polynomial over the original domain.
    pub fn rotate_extended(
        &self,
        poly: &Polynomial<F, ExtendedLagrangeCoeff>,
        rotation: Rotation,
    ) -> Polynomial<F, ExtendedLagrangeCoeff> {
        let new_rotation = ((1 << (self.extended_k - self.k)) * rotation.0.abs()) as usize;

        let mut poly = poly.clone();

        if rotation.0 >= 0 {
            poly.values.rotate_left(new_rotation);
        } else {
            poly.values.rotate_right(new_rotation);
        }

        poly
    }

    /// This takes us from the extended evaluation domain and gets us the
    /// quotient polynomial coefficients.
    ///
    /// This function will panic if the provided vector is not the correct
    /// length.
    // TODO/FIXME: caller should be responsible for truncating
    pub fn extended_to_coeff(&self, mut a: Polynomial<F, ExtendedLagrangeCoeff>) -> Vec<F> {
        assert_eq!(a.values.len(), self.extended_len());

        // Inverse FFT
        Self::ifft(
            &mut a.values,
            &self.evaluation_helper.ext_fft_data,
            true,
            self.extended_ifft_divisor,
        );

        // Distribute powers to move from coset; opposite from the
        // transformation we performed earlier.
        self.distribute_powers_zeta(&mut a.values, false);

        // Truncate it to match the size of the quotient polynomial; the
        // evaluation domain might be slightly larger than necessary because
        // it always lies on a power-of-two boundary.
        a.values
            .truncate((&self.n * self.quotient_poly_degree) as usize);

        a.values
    }

    /// This divides the polynomial (in the extended domain) by the vanishing
    /// polynomial of the $2^k$ size domain.
    pub fn divide_by_vanishing_poly(
        &self,
        mut a: Polynomial<F, ExtendedLagrangeCoeff>,
    ) -> Polynomial<F, ExtendedLagrangeCoeff> {
        assert_eq!(a.values.len(), self.extended_len());

        // Divide to obtain the quotient polynomial in the coset evaluation
        // domain.
        parallelize(&mut a.values, |h, mut index| {
            for h in h {
                h.group_scale(&self.t_evaluations[index % self.t_evaluations.len()]);
                index += 1;
            }
        });

        Polynomial {
            values: a.values,
            _marker: PhantomData,
        }
    }

    /// Given a slice of group elements `[a_0, a_1, a_2, ...]`, this returns
    /// `[a_0, [zeta]a_1, [zeta^2]a_2, a_3, [zeta]a_4, [zeta^2]a_5, a_6, ...]`,
    /// where zeta is a cube root of unity in the multiplicative subgroup with
    /// order (p - 1), i.e. zeta^3 = 1.
    ///
    /// `into_coset` should be set to `true` when moving into the coset,
    /// and `false` when moving out. This toggles the choice of `zeta`.
    fn distribute_powers_zeta(&self, a: &mut [F], into_coset: bool) {
        let coset_powers = if into_coset {
            [self.g_coset, self.g_coset_inv]
        } else {
            [self.g_coset_inv, self.g_coset]
        };
        parallelize(a, |a, mut index| {
            for a in a {
                // Distribute powers to move into/from coset
                let i = index % (coset_powers.len() + 1);
                if i != 0 {
                    a.group_scale(&coset_powers[i - 1]);
                }
                index += 1;
            }
        });
    }

    fn ifft(a: &mut [F], fft_data: &FFTData<F>, is_inv: bool, divisor: F) {
        best_fft(a, &fft_data, is_inv);
        parallelize(a, |a, _| {
            for a in a {
                // Finish iFFT
                a.group_scale(&divisor);
            }
        });
    }

    /// Get the size of the extended domain
    pub fn extended_len(&self) -> usize {
        1 << self.extended_k
    }

    /// Get degree
    pub fn get_degree(&self) -> (u32, u32) {
        (self.k, self.extended_k)
    }

    /// Get $\omega$, the generator of the $2^k$ order multiplicative subgroup.
    pub fn get_omega(&self) -> F {
        self.omega
    }

    /// Get $\omega^{-1}$, the inverse of the generator of the $2^k$ order
    /// multiplicative subgroup.
    pub fn get_omega_inv(&self) -> F {
        self.omega_inv
    }

    /// Get the generator of the extended domain's multiplicative subgroup.
    pub fn get_extended_omega(&self) -> F {
        self.extended_omega
    }

    /// Get the generator of the extended domain's multiplicative subgroup.
    pub fn get_extended_omega_inv(&self) -> F {
        self.extended_omega_inv
    }

    /// Multiplies a value by some power of $\omega$, essentially rotating over
    /// the domain.
    pub fn rotate_omega(&self, value: F, rotation: Rotation) -> F {
        let mut point = value;
        if rotation.0 >= 0 {
            point *= &self.get_omega().pow_vartime(&[rotation.0 as u64]);
        } else {
            point *= &self
                .get_omega_inv()
                .pow_vartime(&[(rotation.0 as i64).abs() as u64]);
        }
        point
    }

    /// Computes evaluations (at the point `x`, where `xn = x^n`) of Lagrange
    /// basis polynomials `l_i(X)` defined such that `l_i(omega^i) = 1` and
    /// `l_i(omega^j) = 0` for all `j != i` at each provided rotation `i`.
    ///
    /// # Implementation
    ///
    /// The polynomial
    ///     $$\prod_{j=0,j \neq i}^{n - 1} (X - \omega^j)$$
    /// has a root at all points in the domain except $\omega^i$, where it evaluates to
    ///     $$\prod_{j=0,j \neq i}^{n - 1} (\omega^i - \omega^j)$$
    /// and so we divide that polynomial by this value to obtain $l_i(X)$. Since
    ///     $$\prod_{j=0,j \neq i}^{n - 1} (X - \omega^j)
    ///       = \frac{X^n - 1}{X - \omega^i}$$
    /// then $l_i(x)$ for some $x$ is evaluated as
    ///     $$\left(\frac{x^n - 1}{x - \omega^i}\right)
    ///       \cdot \left(\frac{1}{\prod_{j=0,j \neq i}^{n - 1} (\omega^i - \omega^j)}\right).$$
    /// We refer to
    ///     $$1 \over \prod_{j=0,j \neq i}^{n - 1} (\omega^i - \omega^j)$$
    /// as the barycentric weight of $\omega^i$.
    ///
    /// We know that for $i = 0$
    ///     $$\frac{1}{\prod_{j=0,j \neq i}^{n - 1} (\omega^i - \omega^j)} = \frac{1}{n}.$$
    ///
    /// If we multiply $(1 / n)$ by $\omega^i$ then we obtain
    ///     $$\frac{1}{\prod_{j=0,j \neq 0}^{n - 1} (\omega^i - \omega^j)}
    ///       = \frac{1}{\prod_{j=0,j \neq i}^{n - 1} (\omega^i - \omega^j)}$$
    /// which is the barycentric weight of $\omega^i$.
    pub fn l_i_range<I: IntoIterator<Item = i32> + Clone>(
        &self,
        x: F,
        xn: F,
        rotations: I,
    ) -> Vec<F> {
        let mut results;
        {
            let rotations = rotations.clone().into_iter();
            results = Vec::with_capacity(rotations.size_hint().1.unwrap_or(0));
            for rotation in rotations {
                let rotation = Rotation(rotation);
                let result = x - self.rotate_omega(F::Scalar::one(), rotation);
                results.push(result);
            }
            results.iter_mut().batch_invert();
        }

        let common = (xn - F::Scalar::one()) * self.barycentric_weight;
        for (rotation, result) in rotations.into_iter().zip(results.iter_mut()) {
            let rotation = Rotation(rotation);
            *result = self.rotate_omega(*result * common, rotation);
        }

        results
    }

    /// Gets the quotient polynomial's degree (as a multiple of n)
    pub fn get_quotient_poly_degree(&self) -> usize {
        self.quotient_poly_degree as usize
    }

    /// Obtain a pinned version of this evaluation domain; a structure with the
    /// minimal parameters needed to determine the rest of the evaluation
    /// domain.
    pub fn pinned(&self) -> PinnedEvaluationDomain<'_, F> {
        PinnedEvaluationDomain {
            k: &self.k,
            extended_k: &self.extended_k,
            omega: &self.omega,
        }
    }
}

/// Represents the minimal parameters that determine an `EvaluationDomain`.
#[allow(dead_code)]
#[derive(Debug)]
pub struct PinnedEvaluationDomain<'a, F: Field> {
    k: &'a u32,
    extended_k: &'a u32,
    omega: &'a F,
}

#[test]
fn test_fft() {
    use crate::poly::{commitment::prev_fft, EvaluationDomain, FFTCache};
    use ark_std::{end_timer, start_timer};
    use pairing::bn256::Fr;
    use rand_core::OsRng;

    for k in 1..20 {
        let rng = OsRng;
        // polynomial degree n = 2^k
        let n = 1u64 << k;
        // polynomial coeffs
        let coeffs: Vec<_> = (0..n).map(|_| Fr::random(rng)).collect();
        // evaluation domain
        let domain: EvaluationDomain<Fr> = EvaluationDomain::new(1, k);
        // fft cache
        let (k, ext_k) = domain.get_degree();
        let fft_cash = FFTHelper::new(
            k as usize,
            ext_k as usize,
            domain.get_omega(),
            domain.get_omega_inv(),
        );

        let mut prev_fft_coeffs = coeffs.clone();
        let mut best_fft_coeffs = coeffs.clone();
        let mut optimized_fft_coeffs = coeffs.clone();

        let message = format!("prev_fft degree {}", k);
        let start = start_timer!(|| message);
        prev_fft(&mut prev_fft_coeffs, domain.get_omega(), k);
        end_timer!(start);

        let message = format!("best_fft degree {}", k);
        let start = start_timer!(|| message);
        best_fft(
            &mut best_fft_coeffs,
            &domain.evaluation_helper.fft_data,
            false,
        );
        end_timer!(start);

        let message = format!("optimized_fft_coeffs degree {}", k);
        let start = start_timer!(|| message);
        fft_cash
            .fft_data
            .bit_reverse
            .sort_bit_reverse(&mut optimized_fft_coeffs);
        fft_cash
            .fft_data
            .butterfly
            .butterfly_arithmetic(&mut optimized_fft_coeffs, n as usize);
        end_timer!(start);

        assert_eq!(prev_fft_coeffs, best_fft_coeffs);
        assert_eq!(optimized_fft_coeffs, best_fft_coeffs);
    }
}

#[test]
fn test_rotate() {
    use rand_core::OsRng;

    use crate::arithmetic::eval_polynomial;
    use pairing::bn256::Fr as Scalar;

    let domain = EvaluationDomain::<Scalar>::new(1, 3);
    let rng = OsRng;

    let mut poly = domain.empty_lagrange();
    assert_eq!(poly.len(), 8);
    for value in poly.iter_mut() {
        *value = Scalar::random(rng);
    }

    let poly_rotated_cur = poly.rotate(Rotation::cur());
    let poly_rotated_next = poly.rotate(Rotation::next());
    let poly_rotated_prev = poly.rotate(Rotation::prev());

    let poly = domain.lagrange_to_coeff(poly);
    let poly_rotated_cur = domain.lagrange_to_coeff(poly_rotated_cur);
    let poly_rotated_next = domain.lagrange_to_coeff(poly_rotated_next);
    let poly_rotated_prev = domain.lagrange_to_coeff(poly_rotated_prev);

    let x = Scalar::random(rng);

    assert_eq!(
        eval_polynomial(&poly[..], x),
        eval_polynomial(&poly_rotated_cur[..], x)
    );
    assert_eq!(
        eval_polynomial(&poly[..], x * domain.omega),
        eval_polynomial(&poly_rotated_next[..], x)
    );
    assert_eq!(
        eval_polynomial(&poly[..], x * domain.omega_inv),
        eval_polynomial(&poly_rotated_prev[..], x)
    );
}

#[test]
fn test_l_i() {
    use rand_core::OsRng;

    use crate::arithmetic::{eval_polynomial, lagrange_interpolate, BaseExt};
    use pairing::bn256::Fr as Scalar;
    let domain = EvaluationDomain::<Scalar>::new(1, 3);

    let mut l = vec![];
    let mut points = vec![];
    for i in 0..8 {
        points.push(domain.omega.pow(&[i, 0, 0, 0]));
    }
    for i in 0..8 {
        let mut l_i = vec![Scalar::zero(); 8];
        l_i[i] = Scalar::one();
        let l_i = lagrange_interpolate(&points[..], &l_i[..]);
        l.push(l_i);
    }

    let x = Scalar::random(OsRng);
    let xn = x.pow(&[8, 0, 0, 0]);

    let evaluations = domain.l_i_range(x, xn, -7..=7);
    for i in 0..8 {
        assert_eq!(eval_polynomial(&l[i][..], x), evaluations[7 + i]);
        assert_eq!(eval_polynomial(&l[(8 - i) % 8][..], x), evaluations[7 - i]);
    }
}
