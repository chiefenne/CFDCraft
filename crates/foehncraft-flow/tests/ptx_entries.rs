//! Every kernel entry the CUDA backend loads (`flow_` + its name, from
//! `foehncraft_kernels::entry::NAMES`) is an entry point of the PTX module
//! (`gpu/ptx-kernels/src/flow.rs`). The module is built apart from the
//! workspace, so a missing line there fails only when a GPU run loads it.

#[test]
fn the_ptx_module_has_every_entry() {
    let path =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../gpu/ptx-kernels/src/flow.rs");
    let src = std::fs::read_to_string(&path).unwrap();
    for name in foehncraft_kernels::entry::NAMES {
        assert!(
            src.contains(&format!("flow_{name} => ")),
            "{} lacks flow_{name}",
            path.display()
        );
    }
}
