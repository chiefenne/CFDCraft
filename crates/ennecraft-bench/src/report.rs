//! G1 report: JSON (machine-readable, for later comparison), Markdown
//! (for docs/benchmarks/g1.md) and an SVG bar chart (for the gallery).

use crate::mesh::Ordering;
use crate::workload::Kernel;
use serde::Serialize;
use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::path::{Path, PathBuf};

/// Complete result of a G1 run.
#[derive(Debug, Clone, Serialize)]
pub struct G1Report {
    /// Version of the code that ran.
    pub ennecraft_version: String,
    /// Date of the run (UTC, YYYY-MM-DD).
    pub date_utc: String,
    /// Operating system and architecture.
    pub platform: String,
    /// Whether the kernel bodies were built with bounds-checked loads.
    pub checked_indexing: bool,
    /// Mesh and matrix sizes.
    pub mesh: MeshInfo,
    /// Timings per backend.
    pub results: Vec<BackendResult>,
    /// G1 verdicts for the candidate/reference pairs that ran.
    pub verdicts: Vec<Verdict>,
}

/// Mesh and matrix sizes.
#[derive(Debug, Clone, Serialize)]
pub struct MeshInfo {
    /// Mesh type (`periodic-box`, `hex`, `tet`, `poly`, `mixed`).
    pub kind: String,
    /// Box dimensions.
    pub dims: (usize, usize, usize),
    /// Cells.
    pub cells: usize,
    /// Faces.
    pub faces: usize,
    /// Cell numbering.
    pub ordering: Ordering,
    /// Non-zeros of the CSR matrix.
    pub csr_nnz: usize,
    /// SELL slice height.
    pub sell_chunk: usize,
    /// SELL sorting window.
    pub sell_sigma: usize,
    /// SELL stored entries including padding.
    pub sell_stored: usize,
}

/// Timings of one backend.
#[derive(Debug, Clone, Serialize)]
pub struct BackendResult {
    /// Backend name.
    pub backend: String,
    /// Device description.
    pub device: String,
    /// Per kernel.
    pub kernels: Vec<KernelResult>,
}

/// Timing and check of one kernel on one backend.
#[derive(Debug, Clone, Serialize)]
pub struct KernelResult {
    /// Kernel.
    pub kernel: Kernel,
    /// Median time, seconds.
    pub median_s: f64,
    /// Minimum time, seconds.
    pub min_s: f64,
    /// Effective bandwidth at the median time, GB/s.
    pub gbytes_per_s: f64,
    /// Bytes of the minimum-traffic model.
    pub bytes: f64,
    /// Largest difference to `cpu-serial`, relative to its largest value.
    pub max_rel_diff: f64,
    /// Bitwise equal to `cpu-serial`.
    pub bitwise: bool,
    /// `max_rel_diff` within the case tolerance.
    pub within_tolerance: bool,
}

/// G1 verdict for one candidate against its reference.
#[derive(Debug, Clone, Serialize)]
pub struct Verdict {
    /// Single-source backend.
    pub candidate: String,
    /// Hand-written reference backend.
    pub reference: String,
    /// Median time of the candidate over the reference, per kernel.
    pub ratios: BTreeMap<Kernel, f64>,
    /// Largest allowed ratio.
    pub limit: f64,
    /// Both backends reproduce the reference results.
    pub results_correct: bool,
    /// Every ratio within the limit and the results correct.
    pub pass: bool,
}

impl G1Report {
    /// Markdown summary.
    pub fn to_markdown(&self) -> String {
        let m = &self.mesh;
        let mut s = String::new();
        let _ = writeln!(s, "## G1 run {} ({})\n", self.date_utc, self.platform);
        let _ = writeln!(
            s,
            "Mesh: {} on {} × {} × {}: {} cells, {} faces, {} ordering. \
             CSR: {} non-zeros. SELL-{}-{}: {} stored entries ({:.1} % padding). \
             Ennecraft {}{}.\n",
            m.kind,
            m.dims.0,
            m.dims.1,
            m.dims.2,
            m.cells,
            m.faces,
            m.ordering.name(),
            m.csr_nnz,
            m.sell_chunk,
            m.sell_sigma,
            m.sell_stored,
            100.0 * (m.sell_stored as f64 / m.csr_nnz as f64 - 1.0),
            self.ennecraft_version,
            if self.checked_indexing {
                ", bounds-checked kernel loads"
            } else {
                ""
            }
        );
        let _ = writeln!(
            s,
            "| Backend | Device | Kernel | Median ms | Min ms | GB/s | Check |"
        );
        let _ = writeln!(s, "|---|---|---|---:|---:|---:|---|");
        for b in &self.results {
            for k in &b.kernels {
                let check = if k.bitwise {
                    "bitwise".to_string()
                } else {
                    format!(
                        "{:.1e}{}",
                        k.max_rel_diff,
                        if k.within_tolerance { "" } else { " FAIL" }
                    )
                };
                let _ = writeln!(
                    s,
                    "| {} | {} | {} | {:.3} | {:.3} | {:.1} | {} |",
                    b.backend,
                    b.device,
                    k.kernel,
                    k.median_s * 1e3,
                    k.min_s * 1e3,
                    k.gbytes_per_s,
                    check
                );
            }
        }
        for v in &self.verdicts {
            let _ = writeln!(
                s,
                "\n**{} against {}: {}** (limit {:.2}; results {})\n",
                v.candidate,
                v.reference,
                if v.pass { "PASS" } else { "FAIL" },
                v.limit,
                if v.results_correct {
                    "correct"
                } else {
                    "NOT correct"
                }
            );
            let _ = writeln!(s, "| Kernel | Time ratio |");
            let _ = writeln!(s, "|---|---:|");
            for (k, q) in &v.ratios {
                let _ = writeln!(s, "| {k} | {q:.3} |");
            }
        }
        if !self.verdicts.iter().any(|v| v.candidate == "cuda-rust") {
            let _ = writeln!(
                s,
                "\nThe GPU pair (cuda-rust against cuda-reference) did not run; \
                 G1 is decided on the GPU."
            );
        }
        s
    }

    /// Bar chart of the effective bandwidth, one group per kernel, one bar
    /// per backend.
    pub fn to_svg(&self) -> String {
        let kernels: Vec<Kernel> = self
            .results
            .first()
            .map(|b| b.kernels.iter().map(|k| k.kernel).collect())
            .unwrap_or_default();
        let nb = self.results.len().max(1);
        let vmax = self
            .results
            .iter()
            .flat_map(|b| b.kernels.iter().map(|k| k.gbytes_per_s))
            .fold(0.0f64, f64::max);
        let top = nice_ceiling(vmax.max(1e-9));
        let (w, h) = (720.0, 400.0);
        let (l, r, t, btm) = (64.0, 16.0, 48.0, 84.0);
        let pw = w - l - r;
        let ph = h - t - btm;
        let group = pw / kernels.len().max(1) as f64;
        let bar = (group * 0.8) / nb as f64;
        // Colour-blind safe categorical palette (Okabe & Ito).
        let colours = [
            "#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00",
        ];
        let mut s = String::new();
        let _ = write!(
            s,
            r##"<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" height="{h}" font-family="Helvetica, Arial, sans-serif" font-size="12">
<rect width="{w}" height="{h}" fill="#ffffff"/>
<text x="{l}" y="22" font-size="15" font-weight="bold" fill="#222">G1 kernel benchmark: effective bandwidth</text>
<text x="{l}" y="38" fill="#555">{} mesh, {} cells, {} ordering, {}, {}</text>
"##,
            self.mesh.kind,
            self.mesh.cells,
            self.mesh.ordering.name(),
            self.platform,
            self.date_utc
        );
        for i in 0..=4 {
            let v = top * i as f64 / 4.0;
            let y = t + ph - ph * v / top;
            let _ = writeln!(
                s,
                r##"<line x1="{l}" x2="{}" y1="{y:.1}" y2="{y:.1}" stroke="#ddd"/><text x="{}" y="{:.1}" text-anchor="end" fill="#555">{}</text>"##,
                l + pw,
                l - 6.0,
                y + 4.0,
                trim(v)
            );
        }
        let _ = writeln!(
            s,
            r##"<text transform="translate(16 {:.1}) rotate(-90)" text-anchor="middle" fill="#333">GB/s</text>"##,
            t + ph / 2.0
        );
        for (gi, k) in kernels.iter().enumerate() {
            let x0 = l + gi as f64 * group + group * 0.1;
            for (bi, b) in self.results.iter().enumerate() {
                if let Some(kr) = b.kernels.iter().find(|x| x.kernel == *k) {
                    let hh = ph * kr.gbytes_per_s / top;
                    let _ = writeln!(
                        s,
                        r##"<rect x="{:.1}" y="{:.1}" width="{:.1}" height="{:.1}" fill="{}"><title>{} {}: {:.1} GB/s</title></rect>"##,
                        x0 + bi as f64 * bar,
                        t + ph - hh,
                        (bar - 2.0).max(1.0),
                        hh,
                        colours[bi % colours.len()],
                        b.backend,
                        k,
                        kr.gbytes_per_s
                    );
                }
            }
            let _ = writeln!(
                s,
                r##"<text x="{:.1}" y="{:.1}" text-anchor="middle" fill="#333">{k}</text>"##,
                l + gi as f64 * group + group / 2.0,
                t + ph + 18.0
            );
        }
        let mut lx = l;
        for (bi, b) in self.results.iter().enumerate() {
            let _ = writeln!(
                s,
                r##"<rect x="{lx:.1}" y="{:.1}" width="12" height="12" fill="{}"/><text x="{:.1}" y="{:.1}" fill="#333">{}</text>"##,
                h - 40.0,
                colours[bi % colours.len()],
                lx + 16.0,
                h - 30.0,
                b.backend
            );
            lx += 24.0 + 7.0 * b.backend.len() as f64;
        }
        s.push_str("</svg>\n");
        s
    }

    /// Writes `NAME.json`, `NAME.md` and `NAME.svg` to `dir`.
    pub fn write(&self, dir: &Path, name: &str) -> anyhow::Result<Vec<PathBuf>> {
        std::fs::create_dir_all(dir)?;
        let files = [
            (
                dir.join(format!("{name}.json")),
                serde_json::to_string_pretty(self)? + "\n",
            ),
            (dir.join(format!("{name}.md")), self.to_markdown()),
            (dir.join(format!("{name}.svg")), self.to_svg()),
        ];
        for (p, text) in &files {
            std::fs::write(p, text)?;
        }
        Ok(files.into_iter().map(|f| f.0).collect())
    }
}

fn nice_ceiling(v: f64) -> f64 {
    let e = 10f64.powf(v.log10().floor());
    for m in [1.0, 2.0, 2.5, 5.0, 10.0] {
        if m * e >= v {
            return m * e;
        }
    }
    10.0 * e
}

fn trim(v: f64) -> String {
    let s = format!("{v:.1}");
    s.strip_suffix(".0").map(str::to_string).unwrap_or(s)
}

/// Today's date in UTC, from the system clock (civil-from-days algorithm of
/// H. Hinnant, *chrono-compatible low-level date algorithms*).
pub fn today_utc() -> String {
    let secs = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_secs() as i64)
        .unwrap_or(0);
    let z = secs.div_euclid(86_400) + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z - era * 146_097;
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = yoe + era * 400 + i64::from(m <= 2);
    format!("{y:04}-{m:02}-{d:02}")
}
