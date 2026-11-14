//! G1 case file (TOML, D39 style: every key optional with a default,
//! unknown keys are errors).

use crate::mesh::Ordering;
use serde::{Deserialize, Serialize};
use std::path::{Path, PathBuf};

/// A G1 benchmark case.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct G1Case {
    /// Synthetic mesh.
    pub mesh: MeshSection,
    /// What to run.
    pub run: RunSection,
    /// SELL-C-σ parameters.
    pub sell: SellSection,
    /// Where to write the report.
    pub output: OutputSection,
}

/// `[mesh]`
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct MeshSection {
    /// Approximate number of cells (a near-cubic periodic box); ignored if
    /// `nx`, `ny` and `nz` are all given.
    pub cells: usize,
    /// Cells in x.
    pub nx: Option<usize>,
    /// Cells in y.
    pub ny: Option<usize>,
    /// Cells in z.
    pub nz: Option<usize>,
    /// `natural` or `shuffled` cell numbering.
    pub ordering: Ordering,
    /// Seed of the random fields and of the shuffled numbering.
    pub seed: u64,
}

impl Default for MeshSection {
    fn default() -> Self {
        Self {
            cells: 10_000_000,
            nx: None,
            ny: None,
            nz: None,
            ordering: Ordering::Natural,
            seed: 2026,
        }
    }
}

impl MeshSection {
    /// Box dimensions.
    pub fn dims(&self) -> (usize, usize, usize) {
        match (self.nx, self.ny, self.nz) {
            (Some(a), Some(b), Some(c)) => (a, b, c),
            _ => crate::mesh::box_dims(self.cells),
        }
    }
}

/// `[run]`
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct RunSection {
    /// Backends, see `cfdcraft bench backends`.
    pub backends: Vec<String>,
    /// Kernels (`face-flux`, `cell-gather`, `spmv-csr`, `spmv-sell`).
    pub kernels: Vec<String>,
    /// Timed repetitions per kernel; the report gives the median and the
    /// minimum.
    pub repeat: usize,
    /// Untimed runs before the timed ones.
    pub warmup: usize,
    /// CPU threads (default: all cores).
    pub threads: Option<usize>,
    /// CUDA device ordinal.
    pub device: usize,
    /// PTX module of the single-source route, built by
    /// `gpu/ptx-kernels/build.sh`.
    pub rust_ptx: PathBuf,
    /// Allowed relative difference between a backend's results and
    /// `cpu-serial` (GPUs may contract a·b+c into one fused operation).
    pub tolerance: f64,
}

impl Default for RunSection {
    fn default() -> Self {
        Self {
            backends: vec![
                "cpu-serial".into(),
                "cpu-threads".into(),
                "cpu-native".into(),
            ],
            kernels: crate::workload::Kernel::ALL
                .iter()
                .map(|k| k.name().to_string())
                .collect(),
            repeat: 10,
            warmup: 2,
            threads: None,
            device: 0,
            rust_ptx: PathBuf::from("target/nvptx64-nvidia-cuda/release/cfdcraft_ptx_kernels.ptx"),
            tolerance: 1e-12,
        }
    }
}

/// `[sell]`
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct SellSection {
    /// Slice height C (32 = one NVIDIA warp).
    pub chunk: usize,
    /// Sorting window σ (1: rows not sorted).
    pub sigma: usize,
}

impl Default for SellSection {
    fn default() -> Self {
        Self {
            chunk: 32,
            sigma: 1,
        }
    }
}

/// `[output]`
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct OutputSection {
    /// Directory of the report files.
    pub directory: PathBuf,
    /// Base name of the report files (`NAME.json`, `NAME.md`, `NAME.svg`).
    pub name: String,
}

impl Default for OutputSection {
    fn default() -> Self {
        Self {
            directory: PathBuf::from("output/g1"),
            name: "g1".into(),
        }
    }
}

impl G1Case {
    /// Reads and validates a case file.
    pub fn from_file(path: &Path) -> anyhow::Result<Self> {
        let text = std::fs::read_to_string(path)
            .map_err(|e| anyhow::anyhow!("cannot read {}: {e}", path.display()))?;
        let case: Self =
            toml::from_str(&text).map_err(|e| anyhow::anyhow!("{}: {e}", path.display()))?;
        case.validate()?;
        Ok(case)
    }

    /// Checks values that the types alone do not.
    pub fn validate(&self) -> anyhow::Result<()> {
        let (a, b, c) = self.mesh.dims();
        anyhow::ensure!(
            a >= 3 && b >= 3 && c >= 3,
            "mesh: at least 3 cells per direction"
        );
        anyhow::ensure!(
            3 * a * b * c < (1usize << 31),
            "mesh: {} cells exceed the 32-bit local index range of one partition",
            a * b * c
        );
        anyhow::ensure!(self.run.repeat >= 1, "run.repeat must be at least 1");
        anyhow::ensure!(!self.run.backends.is_empty(), "run.backends is empty");
        for k in &self.run.kernels {
            anyhow::ensure!(
                crate::workload::Kernel::parse(k).is_some(),
                "run.kernels: unknown kernel `{k}` (face-flux, cell-gather, spmv-csr, spmv-sell)"
            );
        }
        anyhow::ensure!(
            self.sell.chunk >= 1 && self.sell.sigma >= 1,
            "sell.chunk and sell.sigma must be at least 1"
        );
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn defaults_and_unknown_keys() {
        let c: G1Case = toml::from_str("[mesh]\ncells = 1000\n").unwrap();
        assert_eq!(c.run.repeat, 10);
        assert_eq!(c.mesh.dims(), (10, 10, 10));
        c.validate().unwrap();
        assert!(toml::from_str::<G1Case>("[mesh]\ncels = 1000\n").is_err());
        let bad: G1Case = toml::from_str("[run]\nkernels = [\"spmv\"]\n").unwrap();
        assert!(bad.validate().is_err());
    }

    #[test]
    fn shipped_cases_parse() {
        let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/cases/bench");
        let mut n = 0;
        for e in std::fs::read_dir(dir).unwrap() {
            let p = e.unwrap().path();
            if p.extension().is_some_and(|x| x == "toml") {
                G1Case::from_file(&p).unwrap();
                n += 1;
            }
        }
        assert!(n >= 1);
    }
}
