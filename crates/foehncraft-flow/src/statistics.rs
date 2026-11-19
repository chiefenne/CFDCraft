//! Statistics of a transient run (D71): running means, root mean squares
//! and covariances of the fields, means and RMS of the probes and the force
//! coefficients with a confidence interval from batch means (Schmeiser
//! 1982), and phase averages of the probes and coefficients.
//!
//! `[output.statistics]` starts the averaging at `start` (default: the
//! start of the run). Every accepted time step adds its state with the
//! weight of its step Δt, so that variable steps (D65) give time averages.
//! The field moments are updated with the weighted form of Welford's
//! algorithm (West 1979), which keeps the variance free of the
//! cancellation of Σx² − (Σx)²:
//!
//!   W ← W + h,  d = x − m,  m ← m + (h/W) d,  M ← M + h d (x − m),
//!
//! and the covariance of x and y the same with C ← C + h d_x (y − m_y).
//! The mean is m, the variance M/W, the covariance C/W.
//!
//! The series of the probes (u, v, w, p at every probe) and of the force
//! coefficients keep their samples (time, step, value) after `start`. At the
//! end, the averaging window T is cut into `batches` batches of equal time
//! (B, default 20); the time mean of each batch is a batch mean, and the
//! half width of the 95 % confidence interval of the mean is
//! t(0.975, B − 1) s_B/√B with s_B the standard deviation of the batch
//! means. The batch means are only independent when a batch is longer than
//! the correlation time of the series; the von Neumann test of their lag-1
//! correlation tells (Fishman 1978): r₁ = 1 − Σ(b_{i+1} − b_i)²/(2 Σ(b_i −
//! b̄)²) is about normal with mean 0 and variance (B − 2)/(B² − 1) for
//! independent batch means, so r₁ ≤ 1.645 √((B − 2)/(B² − 1)) (one-sided,
//! 95 %; 0.47 for 10 batches, 0.35 for 20) gives the verdict `sufficient`,
//! else `too short` (a longer window, or fewer batches, is needed; STAT001
//! in the report). A negative r₁ (alternating batch means) only widens the
//! interval and passes. A mean known well enough passes too: when the half
//! width widened for the correlation, h √((1 + r₁)/(1 − r₁)), is below
//! `tolerance` (default 10⁻³) of the size √var + |mean| of the series (a
//! slow drift of a well-converged mean).
//!
//! With `phase_period` P the samples are also sorted into `phase_bins`
//! bins of the phase ((t − start)/P mod 1), and the time mean of each bin
//! is written to `statistics_phase.csv` (the phase average of a periodic
//! flow, e.g. a vortex street or a rotor at its blade passing period).

use anyhow::Result;
use serde::{Deserialize, Serialize};

/// `[output.statistics]`.
#[derive(Debug, Clone, Serialize, Deserialize, schemars::JsonSchema, PartialEq)]
#[serde(default, deny_unknown_fields)]
pub struct StatisticsSection {
    /// Simulated time from which the statistics are taken (s); default the
    /// start of the run.
    pub start: Option<f64>,
    /// Running means, RMS and covariances of the fields
    /// (`statistics.vtkhdf`); the probes and coefficients are always
    /// averaged.
    pub fields: bool,
    /// Number of batches of the batch-means confidence interval (≥ 2).
    pub batches: usize,
    /// Period of the phase averages of the probes and coefficients (s);
    /// none without it.
    pub phase_period: Option<f64>,
    /// Number of phase bins per period.
    pub phase_bins: usize,
    /// A mean whose confidence interval, widened for the correlation of
    /// the batch means, is below this fraction of the size of its series
    /// (√var + |mean|) is `sufficient` also when the batch means are
    /// correlated (0: the test of the batch means alone).
    pub tolerance: f64,
}

impl Default for StatisticsSection {
    fn default() -> Self {
        Self {
            start: None,
            fields: true,
            batches: 20,
            phase_period: None,
            phase_bins: 24,
            tolerance: 1e-3,
        }
    }
}

impl StatisticsSection {
    /// Checks the section; `transient` is true with a `[time]` table.
    pub fn validate(&self, transient: bool) -> Result<()> {
        crate::need!(
            transient,
            CASE_INCONSISTENT,
            "[output.statistics] averages over time and needs a [time] table (D71)"
        );
        crate::need!(
            self.batches >= 2,
            CASE_RANGE,
            "[output.statistics] batches must be at least 2 (it is {})",
            self.batches
        );
        crate::need!(
            self.tolerance >= 0.0 && self.tolerance.is_finite(),
            CASE_RANGE,
            "[output.statistics] tolerance must be ≥ 0"
        );
        if let Some(p) = self.phase_period {
            crate::need!(
                p > 0.0 && p.is_finite(),
                CASE_RANGE,
                "[output.statistics] phase_period must be > 0 (s)"
            );
            crate::need!(
                self.phase_bins >= 1,
                CASE_RANGE,
                "[output.statistics] phase_bins must be at least 1"
            );
        }
        if let Some(s) = self.start {
            crate::need!(
                s.is_finite(),
                CASE_RANGE,
                "[output.statistics] start must be finite (s)"
            );
        }
        Ok(())
    }
}

/// Weighted running moments of a set of cell fields (D71).
#[derive(Debug, Clone)]
pub struct FieldMoments {
    /// The names of the variables.
    pub names: Vec<String>,
    /// The pairs (i, j) of the covariances and their names.
    pub pairs: Vec<(usize, usize, String)>,
    /// Total weight W (s).
    pub weight: f64,
    /// Samples.
    pub samples: usize,
    mean: Vec<Vec<f64>>,
    m2: Vec<Vec<f64>>,
    cov: Vec<Vec<f64>>,
    d: Vec<Vec<f64>>,
}

impl FieldMoments {
    /// Moments of `names` on `n` cells with the covariances of `pairs`.
    pub fn new(names: Vec<String>, pairs: Vec<(usize, usize, String)>, n: usize) -> Self {
        let k = names.len();
        Self {
            mean: vec![vec![0.0; n]; k],
            m2: vec![vec![0.0; n]; k],
            cov: vec![vec![0.0; n]; pairs.len()],
            d: vec![vec![0.0; n]; k],
            names,
            pairs,
            weight: 0.0,
            samples: 0,
        }
    }

    /// Adds the values `x` (one slice per variable) with the weight `h`.
    pub fn add(&mut self, h: f64, x: &[&[f64]]) {
        debug_assert_eq!(x.len(), self.names.len());
        if h <= 0.0 {
            return;
        }
        self.weight += h;
        self.samples += 1;
        let r = h / self.weight;
        for (v, xv) in x.iter().enumerate() {
            for (c, &xc) in xv.iter().enumerate() {
                let d = xc - self.mean[v][c];
                self.d[v][c] = d;
                self.mean[v][c] += r * d;
            }
        }
        for (v, xv) in x.iter().enumerate() {
            for (c, &xc) in xv.iter().enumerate() {
                self.m2[v][c] += h * self.d[v][c] * (xc - self.mean[v][c]);
            }
        }
        for (p, (i, j, _)) in self.pairs.iter().enumerate() {
            let (di, mj) = (&self.d[*i], &self.mean[*j]);
            for (c, cv) in self.cov[p].iter_mut().enumerate() {
                *cv += h * di[c] * (x[*j][c] - mj[c]);
            }
        }
    }

    /// The mean of variable `v`.
    pub fn mean(&self, v: usize) -> &[f64] {
        &self.mean[v]
    }

    /// The RMS of the fluctuation of variable `v`, √(M/W).
    pub fn rms(&self, v: usize) -> Vec<f64> {
        let w = self.weight.max(1e-300);
        self.m2[v].iter().map(|m| (m / w).max(0.0).sqrt()).collect()
    }

    /// The covariance of pair `p`, C/W.
    pub fn covariance(&self, p: usize) -> Vec<f64> {
        let w = self.weight.max(1e-300);
        self.cov[p].iter().map(|c| c / w).collect()
    }

    /// The datasets: mean_<name>, rms_<name>, then the covariances.
    pub fn datasets(&self) -> Vec<(String, Vec<f64>)> {
        let mut out = Vec::new();
        for (v, n) in self.names.iter().enumerate() {
            out.push((format!("mean_{n}"), self.mean[v].clone()));
        }
        for (v, n) in self.names.iter().enumerate() {
            out.push((format!("rms_{n}"), self.rms(v)));
        }
        for (p, (_, _, n)) in self.pairs.iter().enumerate() {
            out.push((n.clone(), self.covariance(p)));
        }
        out
    }
}

/// The 0.975 quantile of Student's t distribution with `df` degrees of
/// freedom (table to 30, then the Cornish–Fisher expansion about 1.96).
pub fn t975(df: usize) -> f64 {
    const T: [f64; 30] = [
        12.706, 4.303, 3.182, 2.776, 2.571, 2.447, 2.365, 2.306, 2.262, 2.228, 2.201, 2.179, 2.160,
        2.145, 2.131, 2.120, 2.110, 2.101, 2.093, 2.086, 2.080, 2.074, 2.069, 2.064, 2.060, 2.056,
        2.052, 2.048, 2.045, 2.042,
    ];
    match df {
        0 => f64::INFINITY,
        1..=30 => T[df - 1],
        _ => {
            let z: f64 = 1.959_963_984_540_054;
            let n = df as f64;
            z + (z.powi(3) + z) / (4.0 * n)
                + (5.0 * z.powi(5) + 16.0 * z.powi(3) + 3.0 * z) / (96.0 * n * n)
        }
    }
}

/// The statistics of one series.
#[derive(Debug, Clone, Serialize, PartialEq)]
pub struct SeriesStatistics {
    /// The quantity (`probe.<name>.u`, `coefficients.cd`, …).
    pub quantity: String,
    /// Time mean over the window.
    pub mean: f64,
    /// RMS of the fluctuation about the mean.
    pub rms: f64,
    /// Smallest and largest sample.
    pub min: f64,
    /// See `min`.
    pub max: f64,
    /// Half width of the 95 % confidence interval of the mean from the
    /// batch means.
    pub half_width: f64,
    /// Lag-1 autocorrelation of the batch means.
    pub batch_correlation: f64,
    /// `sufficient` when r₁ passes the von Neumann test (see the module), else
    /// `too short`.
    pub verdict: String,
}

/// The statistics of the series of the probes and coefficients.
#[derive(Debug, Clone, Serialize, PartialEq)]
pub struct StatisticsSummary {
    /// Start and end of the averaging window (s).
    pub window: [f64; 2],
    /// Time steps in the window.
    pub samples: usize,
    /// Batches of the confidence intervals.
    pub batches: usize,
    /// One entry per series.
    pub series: Vec<SeriesStatistics>,
    /// The cell fields of `statistics.vtkhdf`, when written.
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub fields: Vec<String>,
}

/// The samples of the monitored series after the start of the averaging.
#[derive(Debug, Clone, Default)]
pub struct Series {
    /// Names of the series.
    pub names: Vec<String>,
    /// End time of each sample (s).
    pub time: Vec<f64>,
    /// Weight (step) of each sample (s).
    pub step: Vec<f64>,
    /// The values, one row per sample.
    pub values: Vec<Vec<f64>>,
}

impl Series {
    /// Empty series of `names`.
    pub fn new(names: Vec<String>) -> Self {
        Self {
            names,
            ..Default::default()
        }
    }

    /// Adds the sample `v` of the step `h` that ends at `t`.
    pub fn add(&mut self, t: f64, h: f64, v: Vec<f64>) {
        debug_assert_eq!(v.len(), self.names.len());
        if h > 0.0 {
            self.time.push(t);
            self.step.push(h);
            self.values.push(v);
        }
    }

    /// The statistics of every series over the window from `start`, with
    /// `batches` batches of equal time.
    pub fn summary(&self, start: f64, batches: usize) -> StatisticsSummary {
        self.summary_periodic(start, batches, None, 0.0)
    }

    /// As [`Series::summary`]; with a `period` the batches hold whole
    /// periods (as many as fit equally, at least one each, so fewer
    /// batches when the window holds fewer periods), and the rest of the
    /// window enters the mean but no batch: the batch mean of a periodic
    /// signal over part of a period would carry its phase from batch to
    /// batch.
    pub fn summary_periodic(
        &self,
        start: f64,
        batches: usize,
        period: Option<f64>,
        tolerance: f64,
    ) -> StatisticsSummary {
        let end = self.time.last().copied().unwrap_or(start);
        let len = (end - start).max(1e-300);
        let (batches, blen) = match period {
            Some(p) if p > 0.0 && len / p >= 2.0 - 1e-9 => {
                let n = (len / p + 1e-9).floor() as usize;
                let nb = batches.min(n);
                (nb, (n / nb) as f64 * p)
            }
            _ => (batches, len / batches as f64),
        };
        // the size of each series, √(var) + |mean|; the velocity components
        // of a probe share the largest (a component that is round-off in a
        // 2D case then has no correlation to test)
        let size: Vec<f64> = (0..self.names.len())
            .map(|k| {
                let (mut w, mut m, mut m2) = (0.0, 0.0, 0.0);
                for (h, row) in self.step.iter().zip(&self.values) {
                    w += h;
                    let d = row[k] - m;
                    m += h / w * d;
                    m2 += h * d * (row[k] - m);
                }
                (m2 / f64::max(w, 1e-300)).max(0.0).sqrt() + m.abs()
            })
            .collect();
        let group = |k: usize| -> f64 {
            let n = &self.names[k];
            match n.rsplit_once('.') {
                Some((head, "u" | "v" | "w")) => (0..self.names.len())
                    .filter(|&j| {
                        matches!(self.names[j].rsplit_once('.'),
                            Some((h, "u" | "v" | "w")) if h == head)
                    })
                    .map(|j| size[j])
                    .fold(0.0, f64::max),
                _ => size[k],
            }
        };
        let series = (0..self.names.len())
            .map(|k| {
                let (mut w, mut m, mut m2) = (0.0, 0.0, 0.0);
                let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
                let mut bw = vec![0.0; batches];
                let mut bs = vec![0.0; batches];
                for ((t, h), row) in self.time.iter().zip(&self.step).zip(&self.values) {
                    let x = row[k];
                    w += h;
                    let d = x - m;
                    m += h / w * d;
                    m2 += h * d * (x - m);
                    lo = lo.min(x);
                    hi = hi.max(x);
                    // the batch of the middle of the step
                    let mid = t - 0.5 * h;
                    let b = ((mid - start) / blen) as usize;
                    let b = if period.is_none() {
                        b.min(batches - 1)
                    } else {
                        b
                    };
                    if b < batches {
                        bw[b] += h;
                        bs[b] += h * x;
                    }
                }
                let means: Vec<f64> = bw
                    .iter()
                    .zip(&bs)
                    .filter(|(w, _)| **w > 0.0)
                    .map(|(w, s)| s / w)
                    .collect();
                let nb = means.len();
                let scale = group(k);
                let (half, r1) = if nb >= 2 {
                    let mb = means.iter().sum::<f64>() / nb as f64;
                    let var = means.iter().map(|x| (x - mb).powi(2)).sum::<f64>();
                    let sb = (var / (nb - 1) as f64).sqrt();
                    let succ: f64 = means.windows(2).map(|p| (p[1] - p[0]).powi(2)).sum();
                    // batch means equal to round-off: nothing to correlate
                    let r1 = if sb > 1e-12 * scale && var > 0.0 {
                        1.0 - succ / (2.0 * var)
                    } else {
                        0.0
                    };
                    (t975(nb - 1) * sb / (nb as f64).sqrt(), r1)
                } else {
                    (f64::INFINITY, 1.0)
                };
                SeriesStatistics {
                    quantity: self.names[k].clone(),
                    mean: m,
                    rms: (m2 / w.max(1e-300)).max(0.0).sqrt(),
                    min: lo,
                    max: hi,
                    half_width: half,
                    batch_correlation: r1,
                    verdict: if nb >= 3
                        && (r1 <= correlation_bound(nb)
                            || (r1 < 1.0
                                && half * ((1.0 + r1) / (1.0 - r1)).sqrt() <= tolerance * scale))
                    {
                        "sufficient".into()
                    } else {
                        "too short".into()
                    },
                }
            })
            .collect();
        StatisticsSummary {
            window: [start, end],
            samples: self.time.len(),
            batches,
            series,
            fields: Vec::new(),
        }
    }

    /// The phase averages: `bins` bins of the phase ((t − start)/period
    /// mod 1, the middle of each step), the time mean of every series per
    /// bin and the weight in each bin. A bin without samples holds NaN.
    pub fn phase(&self, start: f64, period: f64, bins: usize) -> Vec<(f64, f64, Vec<f64>)> {
        let k = self.names.len();
        let mut w = vec![0.0; bins];
        let mut s = vec![vec![0.0; k]; bins];
        for ((t, h), row) in self.time.iter().zip(&self.step).zip(&self.values) {
            let ph = ((t - 0.5 * h - start) / period).rem_euclid(1.0);
            let b = ((ph * bins as f64) as usize).min(bins - 1);
            w[b] += h;
            for (sv, x) in s[b].iter_mut().zip(row) {
                *sv += h * x;
            }
        }
        (0..bins)
            .map(|b| {
                let centre = (b as f64 + 0.5) / bins as f64;
                let v = if w[b] > 0.0 {
                    s[b].iter().map(|x| x / w[b]).collect()
                } else {
                    vec![f64::NAN; k]
                };
                (centre, w[b], v)
            })
            .collect()
    }
}

/// The end of the statistics of a run: the summary of the series (also in
/// `statistics.csv`, the log and STAT001 for windows that are too short),
/// `statistics_phase.csv` with a phase period, and the field moments in
/// `statistics.vtkhdf`.
#[allow(clippy::too_many_arguments)]
pub(crate) fn finish<F: FnMut(&str)>(
    case: &crate::case::FlowCase,
    cfg: &StatisticsSection,
    start: f64,
    series: &Series,
    fields: Option<&FieldMoments>,
    out: &std::path::Path,
    name: &str,
    mesh: &foehncraft_mesh::Mesh,
    at: (usize, f64),
    lg: &mut crate::case::Log<F>,
    warnings: &mut Vec<crate::diag::Diagnostic>,
) -> Result<StatisticsSummary> {
    let mut sum = series.summary_periodic(start, cfg.batches, cfg.phase_period, cfg.tolerance);
    let pre = |what: &str, units: &str| {
        crate::case::csv_preamble(what, name, Some((at.0, Some(at.1))), units)
    };
    std::fs::write(
        out.join("statistics.csv"),
        summary_csv(
            &pre(
                "statistics",
                "quantity -, mean rms min max (unit of the quantity), half_width_95 (of the \
                 95 % confidence interval of the mean, batch means), batch_correlation - \
                 (lag-1 autocorrelation of the batch means), verdict",
            ),
            &sum,
        ),
    )?;
    if let Some(p) = cfg.phase_period {
        let rows = series.phase(start, p, cfg.phase_bins);
        std::fs::write(
            out.join("statistics_phase.csv"),
            phase_csv(
                &pre(
                    "statistics_phase",
                    "phase - (bin centre, fraction of the period), weight s (time in the \
                     bin), then the time mean of every series in the bin",
                ),
                series,
                &rows,
            ),
        )?;
    }
    lg.line("");
    lg.line(&format!(
        "Statistics (D71) from t = {:.6e} s to {:.6e} s, {} steps, {} batches:",
        sum.window[0], sum.window[1], sum.samples, sum.batches
    ));
    let short: Vec<&str> = sum
        .series
        .iter()
        .filter(|e| e.verdict != "sufficient")
        .map(|e| e.quantity.as_str())
        .collect();
    for e in &sum.series {
        lg.line(&format!(
            "  {:<28} mean {:>13.6e} ± {:.2e}   rms {:.4e}   r1 {:+.2}  {}",
            e.quantity, e.mean, e.half_width, e.rms, e.batch_correlation, e.verdict
        ));
    }
    if !short.is_empty() {
        let m = format!(
            "the averaging window is too short for {} of {} statistics ({}{}): the batch means \
             are correlated",
            short.len(),
            sum.series.len(),
            short[..short.len().min(4)].join(", "),
            if short.len() > 4 { ", …" } else { "" }
        );
        lg.warn(&format!("{}: {m}", crate::diag::STAT_WINDOW));
        warnings.push(crate::diag::Diagnostic::new(crate::diag::STAT_WINDOW, m));
    }
    if let Some(f) = fields {
        let data = f.datasets();
        sum.fields = data.iter().map(|(n, _)| n.clone()).collect();
        if case.output.vtkhdf {
            write_fields(
                &out.join("statistics.vtkhdf"),
                mesh,
                &data,
                name,
                start,
                at,
                f,
            )?;
            lg.note(&format!(
                "statistics of {} fields written (statistics.vtkhdf, {} steps)",
                f.names.len(),
                f.samples
            ));
        }
    }
    Ok(sum)
}

#[cfg(feature = "hdf5")]
fn write_fields(
    file: &std::path::Path,
    mesh: &foehncraft_mesh::Mesh,
    data: &[(String, Vec<f64>)],
    name: &str,
    start: f64,
    at: (usize, f64),
    f: &FieldMoments,
) -> Result<()> {
    use foehncraft_mesh::vtkhdf::{Attr, CellField, write_all};
    let cells: Vec<CellField> = data
        .iter()
        .map(|(n, v)| CellField { name: n, values: v })
        .collect();
    let attrs = [
        ("foehncraft_format", Attr::Str("foehncraft-statistics")),
        ("foehncraft_format_version", Attr::U64(1)),
        ("foehncraft_version", Attr::Str(foehncraft_core::VERSION)),
        ("case", Attr::Str(name)),
        ("step", Attr::U64(at.0 as u64)),
        ("time", Attr::F64(at.1)),
        ("statistics_start", Attr::F64(start)),
        ("statistics_weight", Attr::F64(f.weight)),
        ("statistics_samples", Attr::U64(f.samples as u64)),
    ];
    write_all(file, mesh, &cells, &[], &attrs)
}

#[cfg(not(feature = "hdf5"))]
fn write_fields(
    _: &std::path::Path,
    _: &foehncraft_mesh::Mesh,
    _: &[(String, Vec<f64>)],
    _: &str,
    _: f64,
    _: (usize, f64),
    _: &FieldMoments,
) -> Result<()> {
    Ok(())
}

/// The one-sided 95 % bound of the von Neumann lag-1 correlation of `nb`
/// independent batch means, 1.645 √((nb − 2)/(nb² − 1)).
pub fn correlation_bound(nb: usize) -> f64 {
    let n = nb as f64;
    1.645 * ((n - 2.0).max(0.0) / (n * n - 1.0).max(1.0)).sqrt()
}

/// The CSV text of the phase averages of `series` (D75 preamble included).
pub fn phase_csv(preamble: &str, series: &Series, rows: &[(f64, f64, Vec<f64>)]) -> String {
    use std::fmt::Write as _;
    let mut s = String::from(preamble);
    s.push_str("phase,weight");
    for n in &series.names {
        let _ = write!(s, ",{n}");
    }
    s.push('\n');
    for (ph, w, v) in rows {
        let _ = write!(s, "{ph:e},{w:e}");
        for x in v {
            let _ = write!(s, ",{x:e}");
        }
        s.push('\n');
    }
    s
}

/// The CSV text of the statistics of the series.
pub fn summary_csv(preamble: &str, s: &StatisticsSummary) -> String {
    use std::fmt::Write as _;
    let mut out = String::from(preamble);
    out.push_str("quantity,mean,rms,min,max,half_width_95,batch_correlation,verdict\n");
    for e in &s.series {
        let _ = writeln!(
            out,
            "{},{:e},{:e},{:e},{:e},{:e},{:e},{}",
            e.quantity, e.mean, e.rms, e.min, e.max, e.half_width, e.batch_correlation, e.verdict
        );
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The weighted Welford update gives the time mean, variance and
    /// covariance of samples with unequal weights, to round-off, also for
    /// a large mean (no cancellation).
    #[test]
    fn weighted_moments_match_the_two_pass_values() {
        let n = 200;
        let mut seed = 0x9e37_79b9_7f4a_7c15u64;
        let mut rnd = || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            (seed >> 11) as f64 / (1u64 << 53) as f64
        };
        let h: Vec<f64> = (0..n).map(|_| 0.1 + rnd()).collect();
        let x: Vec<f64> = (0..n).map(|_| 1e6 + rnd()).collect();
        let y: Vec<f64> = x.iter().map(|x| 3.0 * (x - 1e6) + 0.2 * rnd()).collect();
        let mut m = FieldMoments::new(vec!["x".into(), "y".into()], vec![(0, 1, "xy".into())], 1);
        for i in 0..n {
            m.add(h[i], &[&[x[i]], &[y[i]]]);
        }
        let w: f64 = h.iter().sum();
        let mx = h.iter().zip(&x).map(|(h, x)| h * x).sum::<f64>() / w;
        let my = h.iter().zip(&y).map(|(h, y)| h * y).sum::<f64>() / w;
        let vx = h
            .iter()
            .zip(&x)
            .map(|(h, x)| h * (x - mx).powi(2))
            .sum::<f64>()
            / w;
        let cxy = h
            .iter()
            .zip(x.iter().zip(&y))
            .map(|(h, (x, y))| h * (x - mx) * (y - my))
            .sum::<f64>()
            / w;
        assert!((m.mean(0)[0] - mx).abs() < 1e-9 * mx);
        assert!((m.mean(1)[0] - my).abs() < 1e-12);
        assert!((m.rms(0)[0] - vx.sqrt()).abs() < 1e-9 * vx.sqrt());
        assert!((m.covariance(0)[0] - cxy).abs() < 1e-9 * cxy.abs());
        assert_eq!(m.samples, n);
    }

    /// A sine sampled over many periods: mean 0, RMS a/√2; the batches of
    /// whole periods are independent (r₁ small) and the interval holds 0;
    /// a window of a fraction of a period is `too short`.
    #[test]
    fn batch_means_of_a_sine() {
        let a = 2.0;
        let period = 1.0;
        let dt = 0.01;
        let mut s = Series::new(vec!["u".into()]);
        for i in 1..=4000 {
            let t = i as f64 * dt;
            // the time mean of the step, so that the sample is exact
            let v = -a * period / (2.0 * std::f64::consts::PI * dt)
                * ((2.0 * std::f64::consts::PI * t / period).cos()
                    - (2.0 * std::f64::consts::PI * (t - dt) / period).cos());
            s.add(t, dt, vec![v]);
        }
        let r = s.summary(0.0, 20);
        let e = &r.series[0];
        assert!(e.mean.abs() < 1e-12, "{}", e.mean);
        assert!((e.rms - a / 2f64.sqrt()).abs() < 1e-3, "{}", e.rms);
        assert!(e.half_width < 1e-10, "{}", e.half_width);
        assert_eq!(e.verdict, "sufficient");
        // a quarter of a period in 20 batches: strongly correlated
        let mut q = Series::new(vec!["u".into()]);
        for i in 1..=25 {
            let t = i as f64 * dt;
            q.add(t, dt, vec![(2.0 * std::f64::consts::PI * t).sin()]);
        }
        let e = &q.summary(0.0, 20).series[0];
        assert_eq!(e.verdict, "too short", "{e:?}");
        assert!(e.batch_correlation > correlation_bound(20));
        // the phase average of the sine is the sine at the bin centres,
        // to the bin width
        let ph = s.phase(0.0, period, 10);
        for (c, w, v) in &ph {
            assert!(*w > 0.0);
            let exact = a * (2.0 * std::f64::consts::PI * c).sin();
            assert!((v[0] - exact).abs() < 0.07 * a, "{c} {} {exact}", v[0]);
        }
    }

    /// A velocity component that is round-off (w of a 2D case) has no
    /// correlation to test: its verdict follows the size of the probe's
    /// velocity, not its own.
    #[test]
    fn round_off_component_of_a_probe_is_not_too_short() {
        let mut s = Series::new(vec!["probes.a.u".into(), "probes.a.w".into()]);
        for i in 1..=400 {
            let t = i as f64 * 0.01;
            // a slow drift of round-off size in w
            s.add(t, 0.01, vec![1.0 + 0.1 * (40.0 * t).sin(), 1e-19 * t]);
        }
        let e = &s.summary(0.0, 10).series;
        assert_eq!(e[1].verdict, "sufficient", "{e:?}");
        assert_eq!(e[1].batch_correlation, 0.0);
    }

    /// A vortex street with fewer periods than batches would put part of a
    /// period in each batch; with its period the batches hold whole
    /// periods, so the batch means of the sine agree.
    #[test]
    fn batches_hold_whole_periods() {
        let (p, dt) = (0.33, 0.33 / 97.0);
        let mut s = Series::new(vec!["v".into()]);
        let n = (9.3 * p / dt) as usize;
        for i in 1..=n {
            let t = i as f64 * dt;
            s.add(t, dt, vec![(2.0 * std::f64::consts::PI * t / p).sin()]);
        }
        let plain = &s.summary(0.0, 10);
        assert_eq!(plain.series[0].verdict, "too short", "{plain:?}");
        let whole = &s.summary_periodic(0.0, 10, Some(p), 0.0);
        assert_eq!(whole.batches, 9);
        assert_eq!(whole.series[0].verdict, "sufficient", "{whole:?}");
        assert!(whole.series[0].half_width < 1e-3, "{whole:?}");
        // 20 periods in 10 batches: two periods each
        let n = (20.0 * p / dt) as usize + 1;
        let mut s = Series::new(vec!["v".into()]);
        for i in 1..=n {
            let t = i as f64 * dt;
            s.add(t, dt, vec![(2.0 * std::f64::consts::PI * t / p).sin()]);
        }
        let two = &s.summary_periodic(0.0, 10, Some(p), 0.0);
        assert_eq!(two.batches, 10);
        assert_eq!(two.series[0].verdict, "sufficient", "{two:?}");
    }

    /// A converged mean with a slow drift far below its size (the drag of a
    /// vortex street): the batch means are correlated, but the widened
    /// interval is below the tolerance.
    #[test]
    fn tolerance_passes_a_drift_far_below_the_size() {
        let (p, dt) = (0.3, 0.3 / 60.0);
        let mut s = Series::new(vec!["cd".into()]);
        for i in 1..=(20 * 60) {
            let t = i as f64 * dt;
            let x = 3.2 + 0.02 * (4.0 * std::f64::consts::PI * t / p).sin() + 2e-6 * t * t;
            s.add(t, dt, vec![x]);
        }
        let strict = &s.summary_periodic(0.0, 10, Some(p), 0.0).series[0];
        assert_eq!(strict.verdict, "too short", "{strict:?}");
        let tol = &s.summary_periodic(0.0, 10, Some(p), 1e-3).series[0];
        assert_eq!(tol.verdict, "sufficient", "{tol:?}");
    }

    /// White noise: the 95 % interval of the batch means holds the true
    /// mean, and its half width is about 1.96 σ/√N.
    #[test]
    fn batch_means_of_noise() {
        let mut seed = 12345u64;
        let mut rnd = || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            (seed >> 11) as f64 / (1u64 << 53) as f64 - 0.5
        };
        let n = 20_000;
        let mut s = Series::new(vec!["x".into()]);
        for i in 1..=n {
            s.add(i as f64, 1.0, vec![5.0 + rnd()]);
        }
        let e = &s.summary(0.0, 20).series[0];
        let sigma = (1.0f64 / 12.0).sqrt();
        assert!((e.mean - 5.0).abs() <= e.half_width, "{e:?}");
        let expect = 1.96 * sigma / (n as f64).sqrt();
        assert!(
            e.half_width > 0.5 * expect && e.half_width < 2.0 * expect,
            "{e:?}"
        );
        assert_eq!(e.verdict, "sufficient");
        assert!((t975(19) - 2.093).abs() < 1e-12);
        assert!((t975(120) - 1.980).abs() < 2e-3);
    }
}
