//! Binomial checkpointing of a transient run (requirement D59, item 6), for
//! the reverse sweep of a transient adjoint (L17): the adjoint of step i
//! needs the primal state at the start of step i, for i = n − 1, …, 0, but
//! only a few states fit into memory.
//!
//! [`schedule`] gives the order in which a run of n steps with s snapshots
//! (besides the stored start) advances, stores, restores and reverses, with
//! the least number of recomputed steps. It is the binomial scheme of
//! Griewank (1992) and Griewank & Walther (2000): with s snapshots and up to
//! t recomputations of a step, β(s + 1, t) = C(s + 1 + t, t) steps can be
//! reversed. With the smallest t for which β(s + 1, t) ≥ n, the first
//! snapshot of a range of n steps goes to
//! m = max(1, n − β(s, t), β(s + 1, t − 2)): the right part (n − m steps,
//! s − 1 snapshots left) then needs t recomputations at most, the left part
//! (m steps, s snapshots) t − 1, and the total is the least possible,
//! t n − β(s + 2, t − 1) primal steps. The tests compare the count with the
//! optimum of a dynamic program.
//!
//! [`reverse`] runs a schedule on a [`Primal`]; [`Snapshot`] is the primal
//! state of the flow solver in memory (the content of a restart file), so a
//! step recomputed from a snapshot gives the same bits as the first time
//! (the restart of stage 5 is bit-identical).

use crate::simplec::{Flow, State};
use crate::transient::{Clock, OldLevels, Transient};
use anyhow::{Result, ensure};
use std::collections::BTreeMap;

/// One action of a schedule. Level i is the state at the start of step i
/// (level 0 the start of the run, level n its end).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Action {
    /// Run the primal steps `from..to` (the current level is `from`).
    Advance {
        /// Current level.
        from: usize,
        /// Level reached.
        to: usize,
    },
    /// Store the current level.
    Store(usize),
    /// Make a stored level the current one.
    Restore(usize),
    /// Free the snapshot of a level.
    Discard(usize),
    /// Run the adjoint of step i (the current level is i).
    Reverse(usize),
}

/// β(s, t) = C(s + t, t): steps reversible with s − 1 snapshots besides the
/// start and at most t recomputations of each step (saturating).
pub fn beta(s: usize, t: usize) -> usize {
    let mut b: usize = 1;
    for i in 1..=t {
        b = b.saturating_mul(s + i) / i;
    }
    b
}

/// The binomial schedule of `n` steps with `snapshots` snapshots besides
/// the stored start.
pub fn schedule(n: usize, snapshots: usize) -> Vec<Action> {
    let mut a = Vec::new();
    range(0, n, snapshots, &mut a);
    a
}

/// Reverses the steps `l..r`; level `l` is current and stored, `s`
/// snapshots are free.
fn range(l: usize, r: usize, s: usize, a: &mut Vec<Action>) {
    let n = r - l;
    if n == 0 {
        return;
    }
    if n == 1 {
        a.push(Action::Reverse(l));
        return;
    }
    if s == 0 {
        for i in (l..r).rev() {
            if i < r - 1 {
                a.push(Action::Restore(l));
            }
            if i > l {
                a.push(Action::Advance { from: l, to: i });
            }
            a.push(Action::Reverse(i));
        }
        return;
    }
    let mut t = 0;
    while beta(s + 1, t) < n {
        t += 1;
    }
    // optimal positions of the first snapshot: the right part fits into t
    // and needs more than t − 1 recomputations with s − 1 snapshots, the
    // left part fits into t − 1 and needs more than t − 2 with s
    let left_min = if t >= 2 { beta(s + 1, t - 2) } else { 0 };
    let m = l + n.saturating_sub(beta(s, t)).max(left_min).max(1);
    a.push(Action::Advance { from: l, to: m });
    a.push(Action::Store(m));
    range(m, r, s - 1, a);
    a.push(Action::Discard(m));
    a.push(Action::Restore(l));
    range(l, m, s, a);
}

/// What a schedule needs from the primal run.
pub trait Primal {
    /// A stored level.
    type Snapshot;
    /// The current level as a snapshot.
    fn take(&mut self) -> Result<Self::Snapshot>;
    /// Makes a snapshot the current level.
    fn restore(&mut self, s: &Self::Snapshot) -> Result<()>;
    /// Runs step `i` (the current level is `i`).
    fn advance(&mut self, i: usize) -> Result<()>;
}

/// Counts of a reverse sweep.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Stats {
    /// Primal steps run (the first forward sweep included).
    pub steps: usize,
    /// Most snapshots held at once, the start not counted.
    pub max_snapshots: usize,
    /// Adjoint steps.
    pub reversed: usize,
}

/// Runs the schedule of `n` steps with `snapshots` snapshots: stores the
/// current level as the start, then calls `adjoint(p, i)` for
/// i = n − 1, …, 0 with level i current.
pub fn reverse<P: Primal>(
    p: &mut P,
    n: usize,
    snapshots: usize,
    mut adjoint: impl FnMut(&mut P, usize) -> Result<()>,
) -> Result<Stats> {
    let mut store = BTreeMap::new();
    store.insert(0, p.take()?);
    let mut cur = 0;
    let mut st = Stats::default();
    for a in schedule(n, snapshots) {
        match a {
            Action::Advance { from, to } => {
                ensure!(
                    cur == from,
                    "schedule: level {cur} current, {from} expected"
                );
                for i in from..to {
                    p.advance(i)?;
                }
                st.steps += to - from;
                cur = to;
            }
            Action::Store(l) => {
                ensure!(cur == l, "schedule: level {cur} current, {l} stored");
                store.insert(l, p.take()?);
                st.max_snapshots = st.max_snapshots.max(store.len() - 1);
            }
            Action::Restore(l) => {
                let Some(s) = store.get(&l) else {
                    anyhow::bail!("schedule: level {l} is not stored");
                };
                p.restore(s)?;
                cur = l;
            }
            Action::Discard(l) => {
                store.remove(&l);
            }
            Action::Reverse(i) => {
                ensure!(cur == i, "schedule: level {cur} current, step {i} reversed");
                adjoint(p, i)?;
                st.reversed += 1;
            }
        }
    }
    Ok(st)
}

/// The primal state of the flow solver at one level, in host memory: the
/// fields, the old levels and the clock of the time stepping, and the
/// multigrid state (the content of a restart file).
#[derive(Debug, Clone)]
pub struct Snapshot {
    state: State,
    old: OldLevels,
    clock: Clock,
    iteration: usize,
    amg_setups: usize,
    amg_matrix: Option<Vec<f64>>,
}

impl Snapshot {
    /// The current level of `fl` and `tr`.
    pub fn take(fl: &mut Flow, tr: &Transient) -> Result<Self> {
        Ok(Self {
            state: fl.state()?,
            old: tr.old_levels(fl)?,
            clock: tr.clock,
            iteration: fl.iteration(),
            amg_setups: fl.amg_setups,
            amg_matrix: fl.amg_matrix.clone(),
        })
    }

    /// Makes this level the current one of `fl` and `tr`.
    pub fn restore(&self, fl: &mut Flow, tr: &mut Transient) -> Result<()> {
        fl.load(&self.state)?;
        tr.restore(fl, &self.old, self.clock)?;
        fl.iteration = self.iteration;
        match &self.amg_matrix {
            Some(m) => fl.restore_amg(m, self.amg_setups)?,
            None => fl.amg_setups = self.amg_setups,
        }
        Ok(())
    }

    /// Bytes held (for the memory estimate of a schedule).
    pub fn bytes(&self) -> usize {
        let s = &self.state;
        let soa = |v: &foehncraft_fv::disc::Soa3| 8 * (v.x.len() + v.y.len() + v.z.len());
        soa(&s.u)
            + 8 * (s.p.len() + s.mass_flux.len() + s.pb.len())
            + soa(&s.ub)
            + s.grad_u.iter().map(soa).sum::<usize>()
            + self.old.bytes()
            + 8 * self.amg_matrix.as_ref().map_or(0, Vec::len)
    }
}

/// A transient flow run with a fixed time step as a [`Primal`].
pub struct FlowRun<'a> {
    /// Flow solver.
    pub fl: &'a mut Flow,
    /// Time stepping.
    pub tr: &'a mut Transient,
    /// Time step (s).
    pub dt: f64,
}

impl Primal for FlowRun<'_> {
    type Snapshot = Snapshot;
    fn take(&mut self) -> Result<Snapshot> {
        Snapshot::take(self.fl, self.tr)
    }
    fn restore(&mut self, s: &Snapshot) -> Result<()> {
        s.restore(self.fl, self.tr)
    }
    fn advance(&mut self, _i: usize) -> Result<()> {
        self.tr.step(self.fl, self.dt).map(|_| ())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Least primal steps to reverse n steps with s free snapshots (dynamic
    /// program over the position of the first snapshot).
    fn optimum(nmax: usize, smax: usize) -> Vec<Vec<usize>> {
        let mut f = vec![vec![0usize; smax + 1]; nmax + 1];
        for n in 1..=nmax {
            for s in 0..=smax {
                f[n][s] = if n == 1 {
                    0
                } else if s == 0 {
                    n * (n - 1) / 2
                } else {
                    (1..n).map(|m| m + f[n - m][s - 1] + f[m][s]).min().unwrap()
                };
            }
        }
        f
    }

    /// A toy primal: a nonlinear map of a few numbers.
    struct Toy {
        x: Vec<f64>,
        steps: usize,
    }

    impl Primal for Toy {
        type Snapshot = Vec<f64>;
        fn take(&mut self) -> Result<Vec<f64>> {
            Ok(self.x.clone())
        }
        fn restore(&mut self, s: &Vec<f64>) -> Result<()> {
            self.x.clone_from(s);
            Ok(())
        }
        fn advance(&mut self, i: usize) -> Result<()> {
            for (k, v) in self.x.iter_mut().enumerate() {
                *v = (*v * 1.1 + 0.3 * (i + k) as f64).sin();
            }
            self.steps += 1;
            Ok(())
        }
    }

    #[test]
    #[allow(clippy::needless_range_loop)]
    fn schedules_are_optimal_and_reverse_every_level() {
        let (nmax, smax) = (90, 6);
        let opt = optimum(nmax, smax);
        for n in 1..=nmax {
            // the levels of a plain forward run
            let mut toy = Toy {
                x: vec![0.1, 0.2, 0.3],
                steps: 0,
            };
            let mut levels = vec![toy.x.clone()];
            for i in 0..n {
                toy.advance(i).unwrap();
                levels.push(toy.x.clone());
            }
            for s in 0..=smax {
                let mut toy = Toy {
                    x: levels[0].clone(),
                    steps: 0,
                };
                let mut next = n;
                let st = reverse(&mut toy, n, s, |p, i| {
                    assert_eq!(i + 1, next, "steps reversed in order");
                    next = i;
                    assert_eq!(p.x, levels[i], "level {i} of {n} with {s} snapshots");
                    Ok(())
                })
                .unwrap();
                assert_eq!(next, 0);
                assert_eq!(st.reversed, n);
                assert!(st.max_snapshots <= s);
                assert_eq!(
                    st.steps, opt[n][s],
                    "{n} steps, {s} snapshots: {} primal steps, optimum {}",
                    st.steps, opt[n][s]
                );
            }
        }
    }

    #[test]
    fn binomial_counts() {
        // β(s + 1, t) steps with s snapshots and t recomputations
        assert_eq!(beta(1, 1), 2);
        assert_eq!(beta(3, 2), 10);
        assert_eq!(beta(11, 3), 364);
        // 10 000 steps with 20 snapshots: each step runs at most 4 times
        let mut t = 0;
        while beta(21, t) < 10_000 {
            t += 1;
        }
        assert_eq!(t, 4);
        let st = {
            let mut toy = Toy {
                x: vec![0.0],
                steps: 0,
            };
            reverse(&mut toy, 10_000, 20, |_, _| Ok(())).unwrap()
        };
        assert!(st.steps <= 4 * 10_000, "{} steps", st.steps);
    }
}
