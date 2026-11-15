//! Stage-1 acceptance tests on the generated meshes: every cell closed,
//! total volume equal to the domain volume, patches complete.

use inncraft_mesh::Geometry;
use inncraft_mesh::generate::{BoxSpec, Grading, Kind};
use inncraft_mesh::quality::QualityReport;

fn specs() -> Vec<BoxSpec> {
    let mut v = Vec::new();
    for kind in [Kind::Hex, Kind::Tet, Kind::Poly, Kind::Mixed] {
        for dim in [2u8, 3] {
            for perturb in [0.0, 0.25] {
                let mut s = BoxSpec::new(kind, dim, [6, 5, 4]);
                s.lower = [-1.0, 0.0, 0.5];
                s.upper = [2.0, 1.5, 1.0];
                s.perturb = perturb;
                s.grading.y = Grading {
                    ratio: 4.0,
                    symmetric: true,
                };
                s.grading.x = Grading {
                    ratio: 0.5,
                    symmetric: false,
                };
                v.push(s);
            }
        }
    }
    v
}

#[test]
fn generated_meshes_are_closed_and_fill_the_box() {
    for s in specs() {
        let name = format!("{:?} {}D perturb {}", s.kind, s.dimension, s.perturb);
        let m = s.generate().unwrap_or_else(|e| panic!("{name}: {e}"));
        m.validate().unwrap_or_else(|e| panic!("{name}: {e}"));
        let g = Geometry::compute(&m);
        let (q, _) = QualityReport::compute(&m, &g);
        let vol = 3.0 * 1.5 * 0.5;
        assert!(q.ok(), "{name}:\n{}", q.to_markdown(&name));
        assert!(
            ((q.total_volume - vol) / vol).abs() < 1e-12,
            "{name}: volume {}",
            q.total_volume
        );
        assert!(((q.boundary_volume - vol) / vol).abs() < 1e-12, "{name}");
        let area = |p: &str| q.patches.iter().find(|x| x.0 == p).unwrap().2;
        assert!(
            (area("xmin") - 0.75).abs() < 1e-12,
            "{name}: xmin {}",
            area("xmin")
        );
        assert!((area("ymax") - 1.5).abs() < 1e-12, "{name}");
        assert!((area("zmin") - 4.5).abs() < 1e-12, "{name}");
        assert_eq!(q.patches.len(), 6);
        eprintln!(
            "{name}: {} cells, faces/cell {:.0}-{:.0}, shapes {:?}, max non-orth {:.1}",
            q.cells, q.faces_per_cell.min, q.faces_per_cell.max, q.shapes, q.non_orthogonality.max
        );
    }
}
