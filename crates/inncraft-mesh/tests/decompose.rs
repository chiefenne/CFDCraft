//! Partitioning and halo layout: every cell owned once, local meshes valid,
//! volumes add up, the halo exchange delivers the neighbours' values,
//! processor faces match pairwise.

use inncraft_mesh::Geometry;
use inncraft_mesh::decompose::{HaloPlan, decompose};
use inncraft_mesh::generate::{BoxSpec, Kind};
use inncraft_mesh::partition::{Method, PartitionStats, partition};

fn methods() -> Vec<Method> {
    if cfg!(feature = "metis") {
        vec![Method::Rcb, Method::Metis]
    } else {
        vec![Method::Rcb]
    }
}

#[test]
fn decompositions_are_consistent() {
    for kind in [Kind::Hex, Kind::Tet, Kind::Poly, Kind::Mixed] {
        for dim in [2u8, 3] {
            let mut s = BoxSpec::new(kind, dim, [8, 6, 4]);
            s.perturb = 0.2;
            let m = s.generate().unwrap();
            let g = Geometry::compute(&m);
            for method in methods() {
                for parts in [1usize, 2, 3, 5] {
                    let part = partition(&m, &g, parts, method).unwrap();
                    let st = PartitionStats::compute(&m, &part, parts);
                    assert!(st.cells.iter().all(|&c| c > 0));
                    for layers in [1usize, 2] {
                        let name = format!(
                            "{} {dim}D {} {parts} parts {layers} layers",
                            kind.name(),
                            method.name()
                        );
                        let locals = decompose(&m, &part, parts, layers).unwrap();
                        let mut owned = vec![0u32; m.n_cells];
                        let mut vol = 0.0;
                        let mut bfaces = 0;
                        for l in &locals {
                            l.mesh
                                .validate()
                                .unwrap_or_else(|e| panic!("{name}: rank {}: {e}", l.rank));
                            let lg = Geometry::compute(&l.mesh);
                            vol += lg.cell_volume[..l.mesh.n_cells].iter().sum::<f64>();
                            bfaces += l.mesh.n_boundary_faces();
                            for &c in &l.mesh.cell_id[..l.mesh.n_cells] {
                                owned[c as usize] += 1;
                            }
                            // every own cell closed
                            for c in 0..l.mesh.n_cells {
                                let cl = lg.closure(&l.mesh, c);
                                assert!(cl.iter().all(|x| x.abs() < 1e-13), "{name}");
                            }
                            assert_eq!(l.mesh.patches.len(), m.patches.len());
                        }
                        assert!(owned.iter().all(|&k| k == 1), "{name}");
                        assert!(
                            (vol - g.cell_volume.iter().sum::<f64>()).abs() < 1e-12,
                            "{name}"
                        );
                        assert_eq!(bfaces, m.n_boundary_faces(), "{name}");
                        // exchange the global cell index
                        let plans: Vec<&HaloPlan> = locals.iter().map(|l| &l.halo).collect();
                        let mut fields: Vec<Vec<f64>> = locals
                            .iter()
                            .map(|l| {
                                let mut v = vec![-1.0; l.mesh.cell_id.len()];
                                for (x, &id) in v.iter_mut().zip(&l.mesh.cell_id[..l.mesh.n_cells])
                                {
                                    *x = id as f64;
                                }
                                v
                            })
                            .collect();
                        HaloPlan::exchange_all(&plans, &mut fields);
                        for (l, f) in locals.iter().zip(&fields) {
                            for (c, &v) in f.iter().enumerate() {
                                assert_eq!(
                                    v, l.mesh.cell_id[c] as f64,
                                    "{name}: rank {} cell {c}",
                                    l.rank
                                );
                            }
                            // processor faces pair up
                            for (k, nb) in l.halo.neighbours.iter().enumerate() {
                                let other = &locals[nb.rank];
                                let back = other
                                    .halo
                                    .neighbours
                                    .iter()
                                    .position(|x| x.rank == l.rank)
                                    .unwrap();
                                assert_eq!(
                                    l.processor_faces[k], other.processor_faces[back],
                                    "{name}"
                                );
                                assert_eq!(
                                    nb.send.len(),
                                    other.halo.neighbours[back].recv.len(),
                                    "{name}"
                                );
                            }
                        }
                        if layers == 1 {
                            let total: usize = locals
                                .iter()
                                .map(|l| l.processor_faces.iter().sum::<usize>())
                                .sum();
                            assert_eq!(total, 2 * st.cut_faces, "{name}");
                        }
                    }
                }
            }
        }
    }
}
