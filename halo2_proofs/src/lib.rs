//! # halo2_proofs

#![cfg_attr(docsrs, feature(doc_cfg))]
// Build without warnings on stable 1.51 and later.
#![allow(unknown_lints)]
// Disable old lint warnings until our MSRV is at least 1.51.
#![allow(renamed_and_removed_lints)]
// Use the old lint name to build without warnings until our MSRV is at least 1.51.
#![allow(clippy::unknown_clippy_lints)]
// The actual lints we want to disable.
#![allow(
    clippy::op_ref,
    clippy::assign_op_pattern,
    clippy::too_many_arguments,
    clippy::suspicious_arithmetic_impl,
    clippy::many_single_char_names,
    clippy::same_item_push,
    clippy::upper_case_acronyms
)]
#![deny(broken_intra_doc_links)]
#![deny(missing_debug_implementations)]
#![allow(unsafe_code)]
// Remove this once we update pasta_curves
#![allow(unused_imports)]
#![allow(clippy::derive_partial_eq_without_eq)]

pub mod arithmetic;
pub mod circuit;
pub use halo2curves;
mod multicore;
pub mod plonk;
pub mod poly;
pub mod transcript;

pub mod dev;
mod helpers;
pub use helpers::SerdeFormat;

#[cxx::bridge(namespace = "tachyon::halo2")]
mod ffi {
    // Rust types and signatures exposed to C++.
    extern "Rust" {
        type CppG1Affine;
        type CppG1Jacobian;
        type CppFq;
        type CppFr;
    }

    // C++ types and signatures exposed to Rust.
    unsafe extern "C++" {
        include!("halo2_proofs/include/msm.h");

        fn msm(bases: &[CppG1Affine], scalars: &[CppFr]) -> Box<CppG1Jacobian>;
        #[cfg(feature = "tachyon_msm_gpu")]
        fn init_msm_gpu(degree: u8);
        #[cfg(feature = "tachyon_msm_gpu")]
        fn release_msm_gpu();
        #[cfg(feature = "tachyon_msm_gpu")]
        fn msm_gpu(bases: &[CppG1Affine], scalars: &[CppFr]) -> Box<CppG1Jacobian>;
    }
}

#[repr(C)]
#[derive(Debug)]
pub struct CppG1Affine {
    pub x: CppFq,
    pub y: CppFq,
}

#[repr(C)]
#[derive(Debug)]
pub struct CppG1Jacobian {
    pub x: CppFq,
    pub y: CppFq,
    pub z: CppFq,
}

#[repr(transparent)]
#[derive(Debug)]
pub struct CppFq(pub [u64; 4]);

#[repr(transparent)]
#[derive(Debug)]
pub struct CppFr(pub [u64; 4]);
