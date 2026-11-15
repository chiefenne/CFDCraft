//! MPI smoke test, run in CI with several ranks on one machine (O12):
//!
//! ```text
//! cargo build --release -p foehncraft-parallel --features mpi --example mpi_smoke
//! mpirun -n 4 target/release/examples/mpi_smoke
//! ```
//!
//! Checks a ring exchange (the step of a halo exchange) and that the
//! reproducible sum gives the same bits on every rank and in a second call.

use foehncraft_parallel::{Communicator, MpiComm};

fn main() {
    let comm = MpiComm::init().expect("MPI initialisation");
    let (rank, size) = (comm.rank(), comm.size());

    // Ring exchange: send to rank + 1, receive from rank - 1.
    let n = 1000;
    let send: Vec<f64> = (0..n).map(|i| (rank * n + i) as f64).collect();
    let mut recv = vec![-1.0; n];
    let dest = (rank + 1) % size;
    let source = (rank + size - 1) % size;
    comm.send_recv_f64(&send, dest, &mut recv, source);
    let ok_ring = recv
        .iter()
        .enumerate()
        .all(|(i, &v)| v == (source * n + i) as f64);

    // Reproducible sum of values with rounding error.
    let local = (0..100_000)
        .map(|i| foehncraft_core::rng::uniform([i, rank as u64, 0]))
        .sum::<f64>()
        * 1e-3
        + 0.1;
    let s1 = comm.sum_reproducible(local);
    let s2 = comm.sum_reproducible(local);
    let all = comm.all_gather_f64(s1);
    let ok_sum = s1.to_bits() == s2.to_bits() && all.iter().all(|v| v.to_bits() == s1.to_bits());

    let ok = comm.sum_reproducible(if ok_ring && ok_sum { 0.0 } else { 1.0 }) == 0.0;
    if rank == 0 {
        println!(
            "{}",
            MpiComm::library_version().lines().next().unwrap_or("")
        );
        println!("ranks: {size}, ring exchange: {ok_ring}, reproducible sum: {ok_sum} ({s1:.17e})");
        println!(
            "{}",
            if ok {
                "MPI smoke test passed"
            } else {
                "MPI smoke test FAILED"
            }
        );
    }
    drop(comm);
    if !ok {
        std::process::exit(1);
    }
}
