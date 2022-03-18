//! This module provides common utilities, traits and structures for group,
//! field and polynomial arithmetic.

use super::multicore;
use crate::poly::FFTData;
pub use ff::Field;
use group::{
    ff::{BatchInvert, PrimeField},
    Group as _,
};
pub use pairing::arithmetic::*;

fn multiexp_serial<C: CurveAffine>(coeffs: &[C::Scalar], bases: &[C], acc: &mut C::Curve) {
    let coeffs: Vec<_> = coeffs.iter().map(|a| a.to_repr()).collect();

    let c = if bases.len() < 4 {
        1
    } else if bases.len() < 32 {
        3
    } else {
        (f64::from(bases.len() as u32)).ln().ceil() as usize
    };

    fn get_at<F: PrimeField>(segment: usize, c: usize, bytes: &F::Repr) -> usize {
        let skip_bits = segment * c;
        let skip_bytes = skip_bits / 8;

        if skip_bytes >= 32 {
            return 0;
        }

        let mut v = [0; 8];
        for (v, o) in v.iter_mut().zip(bytes.as_ref()[skip_bytes..].iter()) {
            *v = *o;
        }

        let mut tmp = u64::from_le_bytes(v);
        tmp >>= skip_bits - (skip_bytes * 8);
        tmp = tmp % (1 << c);

        tmp as usize
    }

    let segments = (256 / c) + 1;

    for current_segment in (0..segments).rev() {
        for _ in 0..c {
            *acc = acc.double();
        }

        #[derive(Clone, Copy)]
        enum Bucket<C: CurveAffine> {
            None,
            Affine(C),
            Projective(C::Curve),
        }

        impl<C: CurveAffine> Bucket<C> {
            fn add_assign(&mut self, other: &C) {
                *self = match *self {
                    Bucket::None => Bucket::Affine(*other),
                    Bucket::Affine(a) => Bucket::Projective(a + *other),
                    Bucket::Projective(mut a) => {
                        a += *other;
                        Bucket::Projective(a)
                    }
                }
            }

            fn add(self, mut other: C::Curve) -> C::Curve {
                match self {
                    Bucket::None => other,
                    Bucket::Affine(a) => {
                        other += a;
                        other
                    }
                    Bucket::Projective(a) => other + &a,
                }
            }
        }

        let mut buckets: Vec<Bucket<C>> = vec![Bucket::None; (1 << c) - 1];

        for (coeff, base) in coeffs.iter().zip(bases.iter()) {
            let coeff = get_at::<C::Scalar>(current_segment, c, coeff);
            if coeff != 0 {
                buckets[coeff - 1].add_assign(base);
            }
        }

        // Summation by parts
        // e.g. 3a + 2b + 1c = a +
        //                    (a) + b +
        //                    ((a) + b) + c
        let mut running_sum = C::Curve::identity();
        for exp in buckets.into_iter().rev() {
            running_sum = exp.add(running_sum);
            *acc = *acc + &running_sum;
        }
    }
}

/// Performs a small multi-exponentiation operation.
/// Uses the double-and-add algorithm with doublings shared across points.
pub fn small_multiexp<C: CurveAffine>(coeffs: &[C::Scalar], bases: &[C]) -> C::Curve {
    let coeffs: Vec<_> = coeffs.iter().map(|a| a.to_repr()).collect();
    let mut acc = C::Curve::identity();

    // for byte idx
    for byte_idx in (0..32).rev() {
        // for bit idx
        for bit_idx in (0..8).rev() {
            acc = acc.double();
            // for each coeff
            for coeff_idx in 0..coeffs.len() {
                let byte = coeffs[coeff_idx].as_ref()[byte_idx];
                if ((byte >> bit_idx) & 1) != 0 {
                    acc += bases[coeff_idx];
                }
            }
        }
    }

    acc
}

/// Performs a multi-exponentiation operation.
///
/// This function will panic if coeffs and bases have a different length.
///
/// This will use multithreading if beneficial.
pub fn best_multiexp<C: CurveAffine>(coeffs: &[C::Scalar], bases: &[C]) -> C::Curve {
    assert_eq!(coeffs.len(), bases.len());

    let num_threads = multicore::current_num_threads();
    if coeffs.len() > num_threads {
        let chunk = coeffs.len() / num_threads;
        let num_chunks = coeffs.chunks(chunk).len();
        let mut results = vec![C::Curve::identity(); num_chunks];
        multicore::scope(|scope| {
            let chunk = coeffs.len() / num_threads;

            for ((coeffs, bases), acc) in coeffs
                .chunks(chunk)
                .zip(bases.chunks(chunk))
                .zip(results.iter_mut())
            {
                scope.spawn(move |_| {
                    multiexp_serial(coeffs, bases, acc);
                });
            }
        });
        results.iter().fold(C::Curve::identity(), |a, b| a + b)
    } else {
        let mut acc = C::Curve::identity();
        multiexp_serial(coeffs, bases, &mut acc);
        acc
    }
}

/// Performs a radix-$2$ Fast-Fourier Transformation (FFT) on a vector of size
/// $n = 2^k$, when provided `log_n` = $k$ and an element of multiplicative
/// order $n$ called `omega` ($\omega$). The result is that the vector `a`, when
/// interpreted as the coefficients of a polynomial of degree $n - 1$, is
/// transformed into the evaluations of this polynomial at each of the $n$
/// distinct powers of $\omega$. This transformation is invertible by providing
/// $\omega^{-1}$ in place of $\omega$ and dividing each resulting field element
/// by $n$.
///
/// This will use multithreading if beneficial.
pub fn best_fft<G: Group + std::fmt::Debug>(a: &mut [G], omega: G::Scalar, log_n: u32) {
    let threads = multicore::current_num_threads();
    let log_threads = log2_floor(threads);

    if log_n <= log_threads {
        serial_fft(a, omega, log_n);
    } else {
        parallel_fft(a, omega, log_n, log_threads);
    }
}

/// recursive fft
pub fn recursive_fft<F: FieldExt>(input: &mut [F], fft_data: &FFTData<F>) {
    let stash = input.to_vec();
    let mut elements = 32;
    // bit reverse and bottom four layers butterfly arithmetic
    bottom_layers_butterfly_arithmetic(input, &stash, fft_data, fft_data.half / 2);
    for radix in fft_data.stages.iter() {
        let chunk = fft_data.half * 2 / elements;
        for i in 0..chunk / 2 {
            // element and twiddles offset
            let offset = elements / 2;
            let tw_offset = fft_data.half / elements;
            for p in 0..elements / 2 {
                // indexes of this loop
                let first = i * elements * 2 + p;
                let second = first + offset;
                let third = second + offset;
                let fourth = third + offset;

                // twiddle factor for upper side
                let a_tw_idx = tw_offset * 2 * p;
                let second_tw = input[second] * fft_data.f_twiddles[a_tw_idx];
                let fourth_tw = input[fourth] * fft_data.f_twiddles[a_tw_idx];

                // upper side butterfly arithmetic
                let a = input[first] + second_tw;
                let b = input[first] - second_tw;
                let c = input[third] + fourth_tw;
                let d = input[third] - fourth_tw;

                // twiddle factor for bottom side
                let b_tw1_idx = a_tw_idx / 2;
                let b_tw2_idx = b_tw1_idx + fft_data.half / 2;
                let c_tw = c * fft_data.f_twiddles[b_tw1_idx];
                let d_tw = d * fft_data.f_twiddles[b_tw2_idx];

                // bottom side butterfly arithmetic
                input[first] = a + c_tw;
                input[second] = b + d_tw;
                input[third] = a - c_tw;
                input[fourth] = b - d_tw;
            }
        }
        elements *= radix;
    }
}

fn bottom_layers_butterfly_arithmetic<F: FieldExt>(
    input: &mut [F],
    stash: &Vec<F>,
    fft_data: &FFTData<F>,
    div: usize,
) {
    // twiddles factor 16th root of unity
    let tw_offset = fft_data.half / 8;
    let tw_1 = fft_data.f_twiddles[tw_offset];
    let tw_2 = fft_data.f_twiddles[tw_offset * 2];
    let tw_3 = fft_data.f_twiddles[tw_offset * 3];
    let tw_4 = fft_data.f_twiddles[tw_offset * 4];
    let tw_5 = fft_data.f_twiddles[tw_offset * 5];
    let tw_6 = fft_data.f_twiddles[tw_offset * 6];
    let tw_7 = fft_data.f_twiddles[tw_offset * 7];

    for i in 0..div / 4 {
        // decompress bit reverse indexes
        let offset = 16 * i;
        let el_0 = fft_data.indexes[offset];
        let el_1 = fft_data.indexes[offset + 1];
        let el_2 = fft_data.indexes[offset + 2];
        let el_3 = fft_data.indexes[offset + 3];
        let el_4 = fft_data.indexes[offset + 4];
        let el_5 = fft_data.indexes[offset + 5];
        let el_6 = fft_data.indexes[offset + 6];
        let el_7 = fft_data.indexes[offset + 7];
        let el_8 = fft_data.indexes[offset + 8];
        let el_9 = fft_data.indexes[offset + 9];
        let el_10 = fft_data.indexes[offset + 10];
        let el_11 = fft_data.indexes[offset + 11];
        let el_12 = fft_data.indexes[offset + 12];
        let el_13 = fft_data.indexes[offset + 13];
        let el_14 = fft_data.indexes[offset + 14];
        let el_15 = fft_data.indexes[offset + 15];

        // first layer butterfly arithmetic
        let a_a = stash[el_0] + stash[el_1];
        let a_b = stash[el_0] - stash[el_1];
        let a_c = stash[el_2] + stash[el_3];
        let a_d = stash[el_2] - stash[el_3];
        let b_a = stash[el_4] + stash[el_5];
        let b_b = stash[el_4] - stash[el_5];
        let b_c = stash[el_6] + stash[el_7];
        let b_d = stash[el_6] - stash[el_7];
        let c_a = stash[el_8] + stash[el_9];
        let c_b = stash[el_8] - stash[el_9];
        let c_c = stash[el_10] + stash[el_11];
        let c_d = stash[el_10] - stash[el_11];
        let d_a = stash[el_12] + stash[el_13];
        let d_b = stash[el_12] - stash[el_13];
        let d_c = stash[el_14] + stash[el_15];
        let d_d = stash[el_14] - stash[el_15];

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

        // forth layer butterfly airthmetic
        let k_w_1 = k_b * tw_1;
        let k_w_2 = k_c * tw_2;
        let k_w_3 = k_d * tw_3;
        let k_w_4 = l_a * tw_4;
        let k_w_5 = l_b * tw_5;
        let k_w_6 = l_c * tw_6;
        let k_w_7 = l_d * tw_7;

        input[offset] = i_a + k_a;
        input[offset + 1] = i_b + k_w_1;
        input[offset + 2] = i_c + k_w_2;
        input[offset + 3] = i_d + k_w_3;
        input[offset + 4] = j_a + k_w_4;
        input[offset + 5] = j_b + k_w_5;
        input[offset + 6] = j_c + k_w_6;
        input[offset + 7] = j_d + k_w_7;
        input[offset + 8] = i_a - k_a;
        input[offset + 9] = i_b - k_w_1;
        input[offset + 10] = i_c - k_w_2;
        input[offset + 11] = i_d - k_w_3;
        input[offset + 12] = j_a - k_w_4;
        input[offset + 13] = j_b - k_w_5;
        input[offset + 14] = j_c - k_w_6;
        input[offset + 15] = j_d - k_w_7;
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

fn parallel_fft<G: Group + std::fmt::Debug>(
    a: &mut [G],
    omega: G::Scalar,
    log_n: u32,
    log_threads: u32,
) {
    assert!(log_n >= log_threads);

    let num_threads = 1 << log_threads;
    let log_new_n = log_n - log_threads;
    let mut tmp = vec![vec![G::group_zero(); 1 << log_new_n]; num_threads];
    let new_omega = omega.pow_vartime(&[num_threads as u64, 0, 0, 0]);

    multicore::scope(|scope| {
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

/// This evaluates a provided polynomial (in coefficient form) at `point`.
pub fn eval_polynomial<F: Field>(poly: &[F], point: F) -> F {
    // TODO: parallelize?
    poly.iter()
        .rev()
        .fold(F::zero(), |acc, coeff| acc * point + coeff)
}

/// This computes the inner product of two vectors `a` and `b`.
///
/// This function will panic if the two vectors are not the same size.
pub fn compute_inner_product<F: Field>(a: &[F], b: &[F]) -> F {
    // TODO: parallelize?
    assert_eq!(a.len(), b.len());

    let mut acc = F::zero();
    for (a, b) in a.iter().zip(b.iter()) {
        acc += (*a) * (*b);
    }

    acc
}

/// Divides polynomial `a` in `X` by `X - b` with
/// no remainder.
pub fn kate_division<'a, F: Field, I: IntoIterator<Item = &'a F>>(a: I, mut b: F) -> Vec<F>
where
    I::IntoIter: DoubleEndedIterator + ExactSizeIterator,
{
    b = -b;
    let a = a.into_iter();

    let mut q = vec![F::zero(); a.len() - 1];

    let mut tmp = F::zero();
    for (q, r) in q.iter_mut().rev().zip(a.rev()) {
        let mut lead_coeff = *r;
        lead_coeff.sub_assign(&tmp);
        *q = lead_coeff;
        tmp = lead_coeff;
        tmp.mul_assign(&b);
    }

    q
}

/// This simple utility function will parallelize an operation that is to be
/// performed over a mutable slice.
pub fn parallelize<T: Send, F: Fn(&mut [T], usize) + Send + Sync + Clone>(v: &mut [T], f: F) {
    let n = v.len();
    let num_threads = multicore::current_num_threads();
    let mut chunk = (n as usize) / num_threads;
    if chunk < num_threads {
        chunk = n as usize;
    }

    multicore::scope(|scope| {
        for (chunk_num, v) in v.chunks_mut(chunk).enumerate() {
            let f = f.clone();
            scope.spawn(move |_| {
                let start = chunk_num * chunk;
                f(v, start);
            });
        }
    });
}

fn log2_floor(num: usize) -> u32 {
    assert!(num > 0);

    let mut pow = 0;

    while (1 << (pow + 1)) <= num {
        pow += 1;
    }

    pow
}

/// Returns coefficients of an n - 1 degree polynomial given a set of n points
/// and their evaluations. This function will panic if two values in `points`
/// are the same.
pub fn lagrange_interpolate<F: FieldExt>(points: &[F], evals: &[F]) -> Vec<F> {
    assert_eq!(points.len(), evals.len());
    if points.len() == 1 {
        // Constant polynomial
        return vec![evals[0]];
    } else {
        let mut denoms = Vec::with_capacity(points.len());
        for (j, x_j) in points.iter().enumerate() {
            let mut denom = Vec::with_capacity(points.len() - 1);
            for x_k in points
                .iter()
                .enumerate()
                .filter(|&(k, _)| k != j)
                .map(|a| a.1)
            {
                denom.push(*x_j - x_k);
            }
            denoms.push(denom);
        }
        // Compute (x_j - x_k)^(-1) for each j != i
        denoms.iter_mut().flat_map(|v| v.iter_mut()).batch_invert();

        let mut final_poly = vec![F::zero(); points.len()];
        for (j, (denoms, eval)) in denoms.into_iter().zip(evals.iter()).enumerate() {
            let mut tmp: Vec<F> = Vec::with_capacity(points.len());
            let mut product = Vec::with_capacity(points.len() - 1);
            tmp.push(F::one());
            for (x_k, denom) in points
                .iter()
                .enumerate()
                .filter(|&(k, _)| k != j)
                .map(|a| a.1)
                .zip(denoms.into_iter())
            {
                product.resize(tmp.len() + 1, F::zero());
                for ((a, b), product) in tmp
                    .iter()
                    .chain(std::iter::once(&F::zero()))
                    .zip(std::iter::once(&F::zero()).chain(tmp.iter()))
                    .zip(product.iter_mut())
                {
                    *product = *a * (-denom * x_k) + *b * denom;
                }
                std::mem::swap(&mut tmp, &mut product);
            }
            assert_eq!(tmp.len(), points.len());
            assert_eq!(product.len(), points.len() - 1);
            for (final_coeff, interpolation_coeff) in final_poly.iter_mut().zip(tmp.into_iter()) {
                *final_coeff += interpolation_coeff * eval;
            }
        }
        final_poly
    }
}

/// Given roots [a_0, a_1, ... a_n] returns vanishing polynomials
/// (x - a_0) * (x - a_1) * ... * (x - a_n)
pub fn vanishing_polynomial<F: FieldExt>(roots: &[F]) -> Vec<F> {
    fn mul_with<F: FieldExt>(coeffs: Vec<F>, root: &F) -> Vec<F> {
        let mut ret = vec![F::zero(); coeffs.len() + 1];

        for (i, coeff) in coeffs.iter().enumerate() {
            ret[i] -= *coeff * root;
            ret[i + 1] += coeff;
        }

        ret
    }

    let mut coeffs = vec![F::one()];
    for root in roots {
        coeffs = mul_with(coeffs, root);
    }

    coeffs
}

#[cfg(test)]
use rand_core::OsRng;

#[cfg(test)]
use pairing::bn256::Fr as Fp;

#[test]
fn test_lagrange_interpolate() {
    let rng = OsRng;

    let points = (0..5).map(|_| Fp::random(rng)).collect::<Vec<_>>();
    let evals = (0..5).map(|_| Fp::random(rng)).collect::<Vec<_>>();

    for coeffs in 0..5 {
        let points = &points[0..coeffs];
        let evals = &evals[0..coeffs];

        let poly = lagrange_interpolate(points, evals);
        assert_eq!(poly.len(), points.len());

        for (point, eval) in points.iter().zip(evals) {
            assert_eq!(eval_polynomial(&poly, *point), *eval);
        }
    }
}
