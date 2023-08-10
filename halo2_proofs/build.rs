fn main() {
    cxx_build::bridge("src/lib.rs")
        .files([
            "src/msm.cc",
            #[cfg(feature = "tachyon_msm_gpu")]
            "src/msm_gpu.cc",
        ])
        .flag_if_supported("-std=c++17")
        .compile("halo2_proofs");

    let dep_files = vec![
        "src/lib.rs",
        "src/msm.cc",
        #[cfg(feature = "tachyon_msm_gpu")]
        "src/msm_gpu.cc",
        "include/msm.h",
        "include/msm_gpu.h",
    ];
    for file in dep_files {
        println!("cargo:rerun-if-changed={file}");
    }

    println!("cargo:rustc-link-lib=dylib=tachyon");
}
