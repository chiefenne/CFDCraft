//! Small SVG charts of the flow reports: residual histories, convergence
//! over the mesh size, and profiles against reference data.

use std::fmt::Write as _;

/// One series of a chart.
pub struct Series<'a> {
    /// Legend.
    pub label: &'a str,
    /// Points (x, y).
    pub points: &'a [(f64, f64)],
    /// Draw markers instead of a line.
    pub markers: bool,
}

const COLOURS: [&str; 6] = [
    "#2f6fd0", "#d0582f", "#2e9a5a", "#8a4fc8", "#b8901c", "#555555",
];

/// A chart with linear or logarithmic axes.
pub fn chart(
    title: &str,
    xlabel: &str,
    ylabel: &str,
    series: &[Series],
    logx: bool,
    logy: bool,
) -> String {
    let (x0, x1, y0, y1) = (80.0, 600.0, 390.0, 50.0);
    let tx = |v: f64| if logx { v.max(1e-300).log10() } else { v };
    let ty = |v: f64| if logy { v.max(1e-300).log10() } else { v };
    let all = series
        .iter()
        .flat_map(|s| s.points.iter())
        .filter(|p| (!logx || p.0 > 0.0) && (!logy || p.1 > 0.0));
    let (mut xmin, mut xmax, mut ymin, mut ymax) = (
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::INFINITY,
        f64::NEG_INFINITY,
    );
    for &(x, y) in all {
        xmin = xmin.min(tx(x));
        xmax = xmax.max(tx(x));
        ymin = ymin.min(ty(y));
        ymax = ymax.max(ty(y));
    }
    if !xmin.is_finite() {
        (xmin, xmax, ymin, ymax) = (0.0, 1.0, 0.0, 1.0);
    }
    if logx {
        xmin = xmin.floor();
        xmax = xmax.ceil().max(xmin + 1.0);
    }
    if logy {
        ymin = ymin.floor();
        ymax = ymax.ceil().max(ymin + 1.0);
    } else {
        let pad = 0.05 * (ymax - ymin).max(1e-12);
        ymin -= pad;
        ymax += pad;
    }
    if xmax <= xmin {
        xmax = xmin + 1.0;
    }
    let px = |v: f64| x0 + (tx(v) - xmin) / (xmax - xmin) * (x1 - x0);
    let py = |v: f64| y0 - (ty(v) - ymin) / (ymax - ymin) * (y0 - y1);
    let mut s = String::from(
        r##"<svg xmlns="http://www.w3.org/2000/svg" width="640" height="460" viewBox="0 0 640 460" font-family="sans-serif" font-size="12"><rect width="640" height="460" fill="#fff"/>"##,
    );
    let _ = write!(
        s,
        r##"<text x="340" y="28" text-anchor="middle" font-size="15">{title}</text><rect x="{x0}" y="{y1}" width="{}" height="{}" fill="none" stroke="#444"/><text x="340" y="430" text-anchor="middle">{xlabel}</text><text x="18" y="220" text-anchor="middle" transform="rotate(-90 18 220)">{ylabel}</text>"##,
        x1 - x0,
        y0 - y1
    );
    // ticks
    let ticks = |lo: f64, hi: f64, log: bool| -> Vec<(f64, String)> {
        if log {
            let mut v = Vec::new();
            let mut d = lo;
            while d <= hi + 1e-9 {
                v.push((10f64.powf(d), format!("1e{}", d as i32)));
                d += 1.0;
            }
            v
        } else {
            let span = hi - lo;
            let step = 10f64.powf((span / 5.0).log10().floor());
            let step = [1.0, 2.0, 5.0, 10.0]
                .iter()
                .map(|m| m * step)
                .find(|st| span / st <= 6.0)
                .unwrap_or(step * 10.0);
            let decimals = (-step.log10().floor()).max(0.0) as usize;
            let mut v = Vec::new();
            let mut k = (lo / step).ceil();
            while k * step <= hi + 1e-9 * step {
                let t = k * step;
                // integer multiples of the step: no accumulated round-off,
                // and no "-0"
                let l = format!("{:.*}", decimals, t);
                let l = if l.trim_start_matches(['-', '0', '.']).is_empty() {
                    l.trim_start_matches('-').to_string()
                } else {
                    l
                };
                v.push((t, l));
                k += 1.0;
            }
            v
        }
    };
    for (v, l) in ticks(xmin, xmax, logx) {
        let x = px(v);
        let _ = write!(
            s,
            r##"<line x1="{x:.1}" x2="{x:.1}" y1="{y1}" y2="{y0}" stroke="#eee"/><text x="{x:.1}" y="{}" text-anchor="middle">{l}</text>"##,
            y0 + 18.0
        );
    }
    for (v, l) in ticks(ymin, ymax, logy) {
        let y = py(v);
        let _ = write!(
            s,
            r##"<line x1="{x0}" x2="{x1}" y1="{y:.1}" y2="{y:.1}" stroke="#eee"/><text x="{}" y="{:.1}" text-anchor="end">{l}</text>"##,
            x0 - 6.0,
            y + 4.0
        );
    }
    // legend in the corner with the fewest data points under it
    let all_px: Vec<(f64, f64)> = series
        .iter()
        .flat_map(|se| se.points.iter())
        .filter(|p| (!logx || p.0 > 0.0) && (!logy || p.1 > 0.0))
        .map(|&(x, y)| (px(x), py(y)))
        .collect();
    let (lw, lh) = (180.0, 10.0 + 18.0 * series.len() as f64);
    let corners = [
        (x1 - lw, y1),
        (x1 - lw, y0 - lh),
        (x0 + 10.0, y1),
        (x0 + 10.0, y0 - lh),
    ];
    let (lx, ly0) = corners
        .iter()
        .copied()
        .min_by_key(|&(cx, cy)| {
            all_px
                .iter()
                .filter(|&&(x, y)| {
                    x >= cx - 5.0 && x <= cx + lw && y >= cy - 5.0 && y <= cy + lh + 5.0
                })
                .count()
        })
        .unwrap_or(corners[0]);
    for (i, se) in series.iter().enumerate() {
        let c = COLOURS[i % COLOURS.len()];
        let pts: Vec<(f64, f64)> = se
            .points
            .iter()
            .filter(|p| (!logx || p.0 > 0.0) && (!logy || p.1 > 0.0))
            .map(|&(x, y)| (px(x), py(y)))
            .collect();
        if se.markers {
            for (x, y) in &pts {
                let _ = write!(
                    s,
                    r##"<circle cx="{x:.1}" cy="{y:.1}" r="3.5" fill="none" stroke="{c}" stroke-width="1.5"/>"##
                );
            }
        } else if !pts.is_empty() {
            let d: Vec<String> = pts.iter().map(|(x, y)| format!("{x:.1},{y:.1}")).collect();
            let _ = write!(
                s,
                r##"<polyline points="{}" fill="none" stroke="{c}" stroke-width="1.8"/>"##,
                d.join(" ")
            );
        }
        let ly = ly0 + 16.0 + 18.0 * i as f64;
        if se.markers {
            let _ = write!(
                s,
                r##"<circle cx="{}" cy="{ly}" r="3.5" fill="none" stroke="{c}" stroke-width="1.5"/>"##,
                lx + 20.0
            );
        } else {
            let _ = write!(
                s,
                r##"<line x1="{}" x2="{}" y1="{ly}" y2="{ly}" stroke="{c}" stroke-width="2"/>"##,
                lx + 10.0,
                lx + 30.0
            );
        }
        let _ = write!(
            s,
            r##"<text x="{}" y="{}">{}</text>"##,
            lx + 36.0,
            ly + 4.0,
            se.label
        );
    }
    s.push_str("</svg>\n");
    s
}

/// A horizontal bar chart of signed values around zero (the balance at
/// the end of a run): positive bars to the right in blue, negative bars to
/// the left in orange, the values in a column on the right.
pub fn bars(title: &str, xlabel: &str, bars: &[(String, f64)]) -> String {
    let (x0, x1, top, w) = (240.0, 560.0, 50.0, 680.0);
    let row = 26.0;
    let h = top + row * bars.len().max(1) as f64 + 60.0;
    let big = bars
        .iter()
        .map(|b| b.1.abs())
        .filter(|v| v.is_finite())
        .fold(0.0f64, f64::max)
        .max(1e-300);
    let lo = if bars.iter().any(|b| b.1 < 0.0) {
        -big
    } else {
        0.0
    };
    let hi = if bars.iter().any(|b| b.1 > 0.0) {
        big
    } else {
        0.0
    };
    let (lo, hi) = if hi <= lo { (-1.0, 1.0) } else { (lo, hi) };
    let px = |v: f64| x0 + (v - lo) / (hi - lo) * (x1 - x0);
    let y1 = top + row * bars.len().max(1) as f64;
    let mut s = format!(
        r##"<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}" font-family="sans-serif" font-size="12"><rect width="{w}" height="{h}" fill="#fff"/><text x="{}" y="28" text-anchor="middle" font-size="15">{title}</text><text x="{}" y="{}" text-anchor="middle">{xlabel}</text>"##,
        0.5 * w,
        0.5 * (x0 + x1),
        y1 + 40.0
    );
    for (i, (label, v)) in bars.iter().enumerate() {
        // + 0.0: no "-0" for a zero term
        let v = &(*v + 0.0);
        let y = top + row * i as f64;
        let (a, b) = (px(v.min(0.0)), px(v.max(0.0)));
        let c = if *v >= 0.0 { COLOURS[0] } else { COLOURS[1] };
        let _ = write!(
            s,
            r##"<text x="{}" y="{:.1}" text-anchor="end">{label}</text><rect x="{a:.1}" y="{:.1}" width="{:.1}" height="{}" fill="{c}"/><text x="{}" y="{:.1}" text-anchor="end" font-size="11">{v:+.4e}</text>"##,
            x0 - 10.0,
            y + 17.0,
            y + 4.0,
            (b - a).max(0.5),
            row - 8.0,
            w - 12.0,
            y + 17.0
        );
    }
    let z = px(0.0);
    let _ = write!(
        s,
        r##"<line x1="{z:.1}" x2="{z:.1}" y1="{}" y2="{y1}" stroke="#444"/><line x1="{x0}" x2="{x1}" y1="{y1}" y2="{y1}" stroke="#444"/><text x="{x0}" y="{}" text-anchor="middle">{lo:.2e}</text><text x="{x1}" y="{}" text-anchor="middle">{hi:.2e}</text>"##,
        top - 4.0,
        y1 + 16.0,
        y1 + 16.0
    );
    s.push_str("</svg>\n");
    s
}
