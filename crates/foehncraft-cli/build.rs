//! Records the Git commit of the build for the provenance of every run
//! (D85): `FOEHNCRAFT_GIT_COMMIT`, with `-dirty` when the work tree has
//! uncommitted changes, `unknown` outside a Git checkout. Runs `git` at
//! build time only; the solver itself starts no processes (D82).

use std::path::Path;
use std::process::Command;

fn main() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let git = root.join(".git");
    if git.is_dir() {
        for f in ["HEAD", "index", "refs/heads"] {
            println!("cargo:rerun-if-changed={}", git.join(f).display());
        }
    }
    println!("cargo:rerun-if-env-changed=FOEHNCRAFT_GIT_COMMIT");
    let commit = std::env::var("FOEHNCRAFT_GIT_COMMIT")
        .ok()
        .unwrap_or_else(|| {
            let out = |args: &[&str]| {
                Command::new("git")
                    .args(args)
                    .current_dir(&root)
                    .output()
                    .ok()
                    .filter(|o| o.status.success())
                    .map(|o| String::from_utf8_lossy(&o.stdout).trim().to_string())
            };
            match out(&["rev-parse", "HEAD"]) {
                Some(c) if !c.is_empty() => {
                    let dirty = out(&["status", "--porcelain", "--untracked-files=no"])
                        .is_some_and(|s| !s.is_empty());
                    if dirty { format!("{c}-dirty") } else { c }
                }
                _ => "unknown".into(),
            }
        });
    println!("cargo:rustc-env=FOEHNCRAFT_GIT_COMMIT={commit}");
}
