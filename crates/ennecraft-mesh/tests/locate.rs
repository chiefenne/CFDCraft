//! Stage-1 acceptance (D47): point location and cell-to-cell tracking on
//! every mesh type, against a brute-force search.

use ennecraft_core::rng;
use ennecraft_mesh::Geometry;
use ennecraft_mesh::generate::{BoxSpec, Kind};
use ennecraft_mesh::locate::{Locator, TrackEnd};

fn point(seed: u64, i: u64, lo: [f64; 3], hi: [f64; 3]) -> [f64; 3] {
    let mut x = [0.0; 3];
    for d in 0..3 {
        x[d] = lo[d] + (hi[d] - lo[d]) * rng::uniform([seed, i, d as u64]);
    }
    x
}

#[test]
fn locate_and_track_on_all_mesh_types() {
    for kind in [Kind::Hex, Kind::Tet, Kind::Poly, Kind::Mixed] {
        for dim in [2u8, 3] {
            let mut s = BoxSpec::new(kind, dim, [7, 6, 5]);
            s.perturb = 0.25;
            s.upper = [2.0, 1.0, 0.5];
            let m = s.generate().unwrap();
            let g = Geometry::compute(&m);
            let loc = Locator::new(&m, &g);
            let name = format!("{} {}D", kind.name(), dim);
            let n = 2000u64;
            let mut crossed = 0;
            let mut by_walk = 0;
            for i in 0..n {
                let x = point(7, i, s.lower, s.upper);
                // brute force: every point lies in at least one cell
                // (all cells for the first points, the bin search after that)
                let owners: Vec<usize> = if i < 300 {
                    (0..m.n_cells).filter(|&c| loc.contains(c, x)).collect()
                } else {
                    let c = loc.locate_search(x).unwrap();
                    (0..m.n_cells)
                        .filter(|&d| {
                            d == c
                                || (m.cell_face_iter(c).any(|(f, _)| m.across(f, c) == Some(d))
                                    && loc.contains(d, x))
                        })
                        .collect()
                };
                assert!(!owners.is_empty(), "{name}: point {x:?} in no cell");
                assert!(owners.len() <= 2, "{name}: point {x:?} in cells {owners:?}");
                // walk from a random cell
                let hint = (rng::bits([3, i, 0]) % m.n_cells as u64) as usize;
                let t = loc.track(loc.start(hint), x);
                if t.end == TrackEnd::Reached {
                    by_walk += 1;
                    assert!(
                        owners.contains(&t.position.cell),
                        "{name}: walk ends in {} for {x:?}, owners {owners:?}",
                        t.position.cell
                    );
                    crossed += t.faces_crossed;
                }
                let c = loc.locate(x, hint).unwrap();
                assert!(owners.contains(&c));
            }
            assert_eq!(
                by_walk,
                n,
                "{name}: {} walks did not reach the target",
                n - by_walk
            );
            // a particle path: many short steps, then leave the domain
            let mut pos = loc
                .position_in(
                    loc.locate_search([0.5, 0.5, 0.25]).unwrap(),
                    [0.5, 0.5, 0.25],
                )
                .unwrap();
            for i in 0..500u64 {
                let to = point(11, i, [0.0; 3], s.upper);
                let t = loc.track(pos, to);
                assert_eq!(t.end, TrackEnd::Reached, "{name}: step {i}");
                assert!(loc.contains(t.position.cell, to));
                pos = t.position;
            }
            let t = loc.track(pos, [3.0, 0.5, 0.25]);
            match t.end {
                TrackEnd::Boundary { face } => {
                    assert_eq!(m.patches[m.patch_of(face).unwrap()].name, "xmax");
                    assert!((t.position.x[0] - 2.0).abs() < 1e-9, "{:?}", t.position.x);
                }
                e => panic!("{name}: {e:?}"),
            }
            assert!(loc.locate([2.5, 0.5, 0.25], 0).is_none());
            eprintln!(
                "{name}: {} cells, {n} points located by walking, {:.1} faces crossed per walk",
                m.n_cells,
                crossed as f64 / n as f64
            );
        }
    }
}
