//! Build script: with hypre linked as a shared library (feature `hypre`),
//! sets the rpath of the `cfdcraft` binary to hypre's library directory,
//! which the build script of `cfdcraft-fv` reports.

fn main() {
    println!("cargo:rerun-if-env-changed=DEP_CFDCRAFT_HYPRE_LIBDIR");
    if let Ok(dir) = std::env::var("DEP_CFDCRAFT_HYPRE_LIBDIR")
        && std::env::var("CARGO_CFG_WINDOWS").is_err()
    {
        println!("cargo:rustc-link-arg=-Wl,-rpath,{dir}");
    }
}
