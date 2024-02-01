fn main() {
    cxx_build::bridge("src/lib.rs")
        .flag_if_supported("-std=c++17")
        .compile("halo2_proofs");

    let dep_files = vec!["src/lib.rs"];
    for file in dep_files {
        println!("cargo:rerun-if-changed={file}");
    }

    println!("cargo:rustc-link-lib=dylib=tachyon");
}
