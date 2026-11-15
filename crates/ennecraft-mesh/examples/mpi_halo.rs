//! Halo exchange over MPI, run in CI with 2 and 4 ranks on one machine:
//!
//! ```text
//! cargo build --release -p ennecraft-mesh --features mpi --example mpi_halo
//! mpirun -n 4 target/release/examples/mpi_halo
//! ```
//!
//! Every rank generates the same polyhedral mesh, partitions it the same
//! way (deterministic), extracts its part with a two-layer halo and fills
//! its halo cells by the MPI exchange: the halo values must be the global
//! cell numbers, the halo volumes and centroids those of the global mesh,
//! and the reproducible sum of the cell volumes the volume of the box.

use ennecraft_mesh::Geometry;
use ennecraft_mesh::decompose::extract;
use ennecraft_mesh::generate::{BoxSpec, Kind};
use ennecraft_mesh::partition::{Method, partition};
use ennecraft_parallel::{Communicator, MpiComm};

fn main() {
    let comm = MpiComm::init().expect("MPI initialisation");
    let (rank, size) = (comm.rank(), comm.size());
    let mut spec = BoxSpec::new(Kind::Poly, 3, [10, 8, 6]);
    spec.perturb = 0.2;
    let global = spec.generate().expect("mesh");
    let ggeo = Geometry::compute(&global);
    let method = if cfg!(feature = "metis") {
        Method::Metis
    } else {
        Method::Rcb
    };
    let part = partition(&global, &ggeo, size, method).expect("partition");
    let local = extract(&global, &part, rank, 2).expect("local mesh");
    local.mesh.validate().expect("valid local mesh");

    let mut id: Vec<f64> = vec![-1.0; local.mesh.cell_id.len()];
    for (x, &g) in id.iter_mut().zip(&local.mesh.cell_id[..local.mesh.n_cells]) {
        *x = g as f64;
    }
    local.halo.exchange(&comm, &mut id);
    let ok_ids = id
        .iter()
        .zip(&local.mesh.cell_id)
        .all(|(&v, &g)| v == g as f64);

    let geo = local.geometry(&comm);
    let ok_geo = (0..local.mesh.cell_id.len()).all(|c| {
        let g = local.mesh.cell_id[c] as usize;
        (geo.cell_volume[c] - ggeo.cell_volume[g]).abs() < 1e-15
            && (0..3).all(|d| (geo.cell_centroid[c][d] - ggeo.cell_centroid[g][d]).abs() < 1e-14)
    });
    let vol = comm.sum_reproducible(ennecraft_core::repro::pairwise_sum(
        &geo.cell_volume[..local.mesh.n_cells],
    ));
    let ok_vol = (vol - 1.0).abs() < 1e-12;

    println!(
        "rank {rank}/{size}: {} cells, {} halo cells, neighbours {:?}, ids {}, geometry {}, volume {vol:.15}",
        local.mesh.n_cells,
        local.mesh.n_halo,
        local
            .halo
            .neighbours
            .iter()
            .map(|n| n.rank)
            .collect::<Vec<_>>(),
        ok_ids,
        ok_geo
    );
    let ok = comm.sum_reproducible(if ok_ids && ok_geo && ok_vol { 0.0 } else { 1.0 }) == 0.0;
    if rank == 0 {
        println!(
            "{}",
            if ok {
                "MPI halo exchange: PASS"
            } else {
                "MPI halo exchange: FAIL"
            }
        );
    }
    if !ok {
        std::process::exit(1);
    }
}
