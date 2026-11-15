//! Build script: with the feature `hypre`, compiles the C shim over hypre's
//! public API (`hypre/shim.c`) and links hypre from `HYPRE_DIR`, the install
//! prefix of a sequential hypre build (cmake -DHYPRE_ENABLE_MPI=OFF). A
//! shared library (`libHYPRE.so`, `.dylib`) is preferred and found at run
//! time through an rpath (dependents such as `inncraft-cli` set it from
//! `DEP_INNCRAFT_HYPRE_LIBDIR`); a GPU build of hypre must be shared.

fn main() {
    println!("cargo:rerun-if-changed=hypre/shim.c");
    println!("cargo:rerun-if-env-changed=HYPRE_DIR");
    #[cfg(feature = "hypre")]
    hypre();
}

#[cfg(feature = "hypre")]
fn hypre() {
    use std::path::PathBuf;
    let dir = PathBuf::from(std::env::var("HYPRE_DIR").unwrap_or_else(|_| {
        panic!(
            "the feature `hypre` needs HYPRE_DIR, the install prefix of a sequential \
             hypre build (see docs/recipes.md, recipe 17)"
        )
    }));
    let inc = dir.join("include");
    assert!(
        inc.join("HYPRE.h").exists(),
        "HYPRE_DIR={} has no include/HYPRE.h",
        dir.display()
    );
    let lib = ["lib", "lib64"]
        .iter()
        .map(|l| dir.join(l))
        .find(|l| l.exists())
        .unwrap_or_else(|| panic!("HYPRE_DIR={} has no lib directory", dir.display()));
    cc::Build::new()
        .file("hypre/shim.c")
        .include(&inc)
        .flag_if_supported("-std=c11")
        .warnings(true)
        .compile("inncraft_hypre_shim");
    println!("cargo:rustc-link-search=native={}", lib.display());
    let shared = ["libHYPRE.so", "libHYPRE.dylib"]
        .iter()
        .any(|f| lib.join(f).exists());
    if shared {
        println!("cargo:rustc-link-lib=dylib=HYPRE");
        // rpath for the tests of this crate; dependents read `libdir`
        if !cfg!(windows) {
            println!("cargo:rustc-link-arg=-Wl,-rpath,{}", lib.display());
        }
        println!("cargo:libdir={}", lib.display());
    } else {
        println!("cargo:rustc-link-lib=static=HYPRE");
        if !cfg!(windows) {
            println!("cargo:rustc-link-lib=m");
        }
    }
}
