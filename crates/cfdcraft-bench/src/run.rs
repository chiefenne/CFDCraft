//! The G1 run: build the workload, time every kernel on every backend,
//! check the results, compute the verdict.

use crate::backend::{self, BackendOptions};
use crate::case::G1Case;
use crate::report::{BackendResult, G1Report, KernelResult, MeshInfo, Verdict};
use crate::workload::{Kernel, Workload};
use std::collections::BTreeMap;

/// G1 threshold: the single-source route may be at most this much slower
/// than the hand-written reference (D8).
pub const G1_LIMIT: f64 = 1.10;

/// Candidate and reference backend pairs judged by G1.
pub const PAIRS: [(&str, &str); 2] = [
    ("cpu-threads", "cpu-native"),
    ("cuda-rust", "cuda-reference"),
];

/// Runs a G1 case. `log` receives progress lines.
pub fn run_g1(case: &G1Case, log: &mut dyn FnMut(&str)) -> anyhow::Result<G1Report> {
    case.validate()?;
    let dims = case.mesh.dims();
    log(&format!(
        "building workload: {} mesh on {} x {} x {}, {} ordering",
        case.mesh.kind.name(),
        dims.0,
        dims.1,
        dims.2,
        case.mesh.ordering.name()
    ));
    let w = match case.mesh.kind.generator() {
        None => Workload::new(
            dims,
            case.mesh.ordering,
            case.mesh.seed,
            case.sell.chunk,
            case.sell.sigma,
        ),
        Some(kind) => Workload::stage1(
            kind,
            dims,
            case.mesh.ordering,
            case.mesh.seed,
            case.sell.chunk,
            case.sell.sigma,
        )?,
    };
    log(&format!(
        "{} cells, {} interior faces, {:.2} faces per cell",
        w.mesh.n_cells,
        w.mesh.n_faces(),
        2.0 * w.mesh.n_faces() as f64 / w.mesh.n_cells as f64
    ));
    let kernels: Vec<Kernel> = case
        .run
        .kernels
        .iter()
        .filter_map(|k| Kernel::parse(k))
        .collect();

    // Reference results: the kernel bodies in a sequential loop.
    let mut reference: BTreeMap<Kernel, Vec<f64>> = BTreeMap::new();
    {
        let mut b = backend::create("cpu-serial", &BackendOptions::default())?;
        b.setup(&w)?;
        for &k in &kernels {
            b.run(&w, k)?;
            reference.insert(k, b.output(&w, k)?);
        }
    }

    let opts = BackendOptions {
        threads: case.run.threads,
        device: case.run.device,
        rust_ptx: case.run.rust_ptx.clone(),
    };
    let mut results = Vec::new();
    for name in &case.run.backends {
        let mut b = backend::create(name, &opts)?;
        log(&format!("backend {name}: {}", b.device()));
        b.setup(&w)?;
        let mut kr = Vec::new();
        for &k in &kernels {
            for _ in 0..case.run.warmup {
                b.run(&w, k)?;
            }
            let mut t: Vec<f64> = (0..case.run.repeat)
                .map(|_| b.run(&w, k))
                .collect::<anyhow::Result<_>>()?;
            t.sort_by(f64::total_cmp);
            let median = if t.len() % 2 == 1 {
                t[t.len() / 2]
            } else {
                0.5 * (t[t.len() / 2 - 1] + t[t.len() / 2])
            };
            let out = b.output(&w, k)?;
            let (max_rel_diff, bitwise) = compare(&out, &reference[&k]);
            let bytes = w.bytes(k);
            let r = KernelResult {
                kernel: k,
                median_s: median,
                min_s: t[0],
                gbytes_per_s: bytes / median * 1e-9,
                bytes,
                max_rel_diff,
                bitwise,
                within_tolerance: max_rel_diff <= case.run.tolerance,
            };
            log(&format!(
                "  {:<12} median {:>9.3} ms  min {:>9.3} ms  {:>7.1} GB/s  {}",
                k.name(),
                r.median_s * 1e3,
                r.min_s * 1e3,
                r.gbytes_per_s,
                if bitwise {
                    "bitwise equal".to_string()
                } else {
                    format!("max rel. diff {max_rel_diff:.1e}")
                }
            ));
            kr.push(r);
        }
        results.push(BackendResult {
            backend: name.clone(),
            device: b.device(),
            kernels: kr,
        });
    }

    let verdicts = PAIRS
        .iter()
        .filter_map(|&(cand, refb)| {
            let c = results.iter().find(|r| r.backend == cand)?;
            let r = results.iter().find(|r| r.backend == refb)?;
            let ratios: BTreeMap<Kernel, f64> = c
                .kernels
                .iter()
                .filter_map(|ck| {
                    let rk = r.kernels.iter().find(|x| x.kernel == ck.kernel)?;
                    Some((ck.kernel, ck.median_s / rk.median_s))
                })
                .collect();
            let correct = c
                .kernels
                .iter()
                .chain(&r.kernels)
                .all(|k| k.within_tolerance);
            let pass = correct && ratios.values().all(|&q| q <= G1_LIMIT);
            Some(Verdict {
                candidate: cand.into(),
                reference: refb.into(),
                ratios,
                limit: G1_LIMIT,
                results_correct: correct,
                pass,
            })
        })
        .collect();

    Ok(G1Report {
        cfdcraft_version: cfdcraft_core::VERSION.into(),
        date_utc: crate::report::today_utc(),
        platform: crate::sysinfo::platform(),
        checked_indexing: cfg!(feature = "checked-indexing"),
        mesh: MeshInfo {
            kind: case.mesh.kind.name().to_string(),
            dims,
            cells: w.mesh.n_cells,
            faces: w.mesh.n_faces(),
            ordering: w.ordering,
            csr_nnz: w.csr.nnz(),
            sell_chunk: w.sell.chunk,
            sell_sigma: w.sell.sigma,
            sell_stored: w.sell.stored(),
        },
        results,
        verdicts,
    })
}

/// Largest difference relative to the largest reference value, and whether
/// the two are bitwise equal.
fn compare(out: &[f64], reference: &[f64]) -> (f64, bool) {
    if out.len() != reference.len() {
        return (f64::INFINITY, false);
    }
    let scale = reference
        .iter()
        .fold(0.0f64, |m, v| m.max(v.abs()))
        .max(f64::MIN_POSITIVE);
    let mut d = 0.0f64;
    let mut bitwise = true;
    for (a, b) in out.iter().zip(reference) {
        bitwise &= a.to_bits() == b.to_bits();
        d = d.max((a - b).abs());
    }
    (d / scale, bitwise)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn small_cpu_run_is_bitwise_reproducible() {
        let mut case = G1Case::default();
        case.mesh.cells = 20_000;
        case.mesh.ordering = crate::mesh::Ordering::Shuffled;
        case.sell.sigma = 8;
        case.run.repeat = 2;
        case.run.warmup = 0;
        case.run.threads = Some(3);
        let rep = run_g1(&case, &mut |_| {}).unwrap();
        assert_eq!(rep.results.len(), 3);
        for b in &rep.results {
            for k in &b.kernels {
                assert!(k.bitwise, "{} {}", b.backend, k.kernel);
            }
        }
        assert_eq!(rep.verdicts.len(), 1);
        assert!(rep.verdicts[0].results_correct);
    }
}
