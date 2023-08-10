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
pub mod ffi {
    // Rust types and signatures exposed to C++.
    extern "Rust" {
        type G1MSM;
        type G1MSMGpu;
        type G1Point2;
        type G1JacobianPoint;
        type Fq;
        type Fr;
    }

    // C++ types and signatures exposed to Rust.
    unsafe extern "C++" {
        include!("halo2_proofs/include/msm.h");
        include!("halo2_proofs/include/msm_gpu.h");

        fn create_g1_msm(degree: u8) -> Box<G1MSM>;
        fn destroy_g1_msm(msm: Box<G1MSM>);
        unsafe fn g1_msm(
            msm: *mut G1MSM,
            bases: &[G1Point2],
            scalars: &[Fr],
        ) -> Box<G1JacobianPoint>;
        #[cfg(feature = "tachyon_msm_gpu")]
        fn create_g1_msm_gpu(degree: u8, algorithm: i32) -> Box<G1MSMGpu>;
        #[cfg(feature = "tachyon_msm_gpu")]
        fn destroy_g1_msm_gpu(msm: Box<G1MSMGpu>);
        #[cfg(feature = "tachyon_msm_gpu")]
        unsafe fn g1_msm_gpu(
            msm: *mut G1MSMGpu,
            bases: &[G1Point2],
            scalars: &[Fr],
        ) -> Box<G1JacobianPoint>;
    }
}

#[derive(Debug)]
pub struct G1MSM;

#[derive(Debug)]
pub struct G1MSMGpu;

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

lazy_static::lazy_static! {
    pub static ref MSM_GPU: std::sync::Mutex<stdint::uintptr_t> =
        std::sync::Mutex::new(0);
}
