fn main() {
    let src_files = [
        "src/bn254_blake2b_writer.cc",
        "src/bn254_evals.cc",
        "src/bn254_poly.cc",
        "src/bn254_proving_key.cc",
        "src/bn254_rational_evals.cc",
        "src/bn254_shplonk_prover.cc",
        "src/xor_shift_rng.cc",
    ];
    cxx_build::bridges(["src/bn254.rs", "src/xor_shift_rng.rs"])
        .files(src_files)
        .flag_if_supported("-std=c++17")
        .compile("halo2_proofs");

    let mut dep_files = vec![
        "include/bn254_blake2b_writer.h",
        "include/bn254_evals.h",
        "include/bn254_poly.h",
        "include/bn254_proving_key.h",
        "include/bn254_rational_evals.h",
        "include/bn254_shplonk_prover.h",
        "include/xor_shift_rng.h",
        "src/bn254.rs",
        "src/rust_vec.h",
        "src/xor_shift_rng.rs",
    ];
    dep_files.extend_from_slice(&src_files);
    for file in dep_files {
        println!("cargo:rerun-if-changed={file}");
    }

    println!("cargo:rustc-link-lib=dylib=tachyon");
}
