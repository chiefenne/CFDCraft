//! Communicator abstraction.

use ennecraft_core::repro::pairwise_sum;

/// Communication between the ranks of a run.
///
/// Reductions are reproducible (D31): every rank receives the same bits, and
/// a rerun on the same number of ranks gives the same bits. The value can
/// change with the number of ranks, because the summation order follows the
/// decomposition.
pub trait Communicator {
    /// Rank of this process, `0..size()`.
    fn rank(&self) -> usize;
    /// Number of ranks.
    fn size(&self) -> usize;
    /// Waits until every rank has arrived.
    fn barrier(&self);
    /// One value from every rank, in rank order, on every rank.
    fn all_gather_f64(&self, local: f64) -> Vec<f64>;
    /// Sends `send` to rank `dest` and receives a message of the same length
    /// from rank `source` into `recv`. The basic step of a halo exchange.
    fn send_recv_f64(&self, send: &[f64], dest: usize, recv: &mut [f64], source: usize);

    /// Sum over all ranks, bitwise identical on every rank: the local values
    /// are gathered in rank order and summed in the fixed pairwise order.
    fn sum_reproducible(&self, local: f64) -> f64 {
        pairwise_sum(&self.all_gather_f64(local))
    }

    /// Maximum over all ranks.
    fn max(&self, local: f64) -> f64 {
        self.all_gather_f64(local)
            .into_iter()
            .fold(f64::NEG_INFINITY, f64::max)
    }

    /// Name of the implementation, for logs.
    fn name(&self) -> &'static str;
}

/// The communicator of a run on one rank.
#[derive(Debug, Default, Clone, Copy)]
pub struct SerialComm;

impl Communicator for SerialComm {
    fn rank(&self) -> usize {
        0
    }
    fn size(&self) -> usize {
        1
    }
    fn barrier(&self) {}
    fn all_gather_f64(&self, local: f64) -> Vec<f64> {
        vec![local]
    }
    fn send_recv_f64(&self, send: &[f64], dest: usize, recv: &mut [f64], source: usize) {
        assert_eq!((dest, source), (0, 0), "a serial run has only rank 0");
        recv.copy_from_slice(send);
    }
    fn name(&self) -> &'static str {
        "serial"
    }
}

#[cfg(feature = "mpi")]
pub mod mpi_comm {
    //! MPI implementation of [`Communicator`].

    use super::Communicator;
    use mpi::environment::Universe;
    use mpi::topology::{Communicator as _, SimpleCommunicator};
    use mpi::traits::*;

    /// Communicator over `MPI_COMM_WORLD`. MPI is finalised when the value
    /// is dropped.
    pub struct MpiComm {
        _universe: Universe,
        world: SimpleCommunicator,
    }

    impl MpiComm {
        /// Initialises MPI. Returns `None` if MPI was already initialised.
        pub fn init() -> Option<Self> {
            let universe = mpi::initialize()?;
            let world = universe.world();
            Some(Self {
                _universe: universe,
                world,
            })
        }

        /// Name and version of the MPI library.
        pub fn library_version() -> String {
            mpi::environment::library_version().unwrap_or_default()
        }
    }

    impl Communicator for MpiComm {
        fn rank(&self) -> usize {
            self.world.rank() as usize
        }
        fn size(&self) -> usize {
            self.world.size() as usize
        }
        fn barrier(&self) {
            self.world.barrier();
        }
        fn all_gather_f64(&self, local: f64) -> Vec<f64> {
            let mut out = vec![0.0; self.size()];
            self.world.all_gather_into(&local, &mut out[..]);
            out
        }
        fn send_recv_f64(&self, send: &[f64], dest: usize, recv: &mut [f64], source: usize) {
            let d = self.world.process_at_rank(dest as i32);
            let s = self.world.process_at_rank(source as i32);
            mpi::request::scope(|scope| {
                let sreq = d.immediate_send(scope, send);
                s.receive_into(recv);
                sreq.wait();
            });
        }
        fn name(&self) -> &'static str {
            "mpi"
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn serial_reductions() {
        let c = SerialComm;
        assert_eq!(c.sum_reproducible(2.5), 2.5);
        assert_eq!(c.max(-1.0), -1.0);
        let mut r = [0.0; 3];
        c.send_recv_f64(&[1.0, 2.0, 3.0], 0, &mut r, 0);
        assert_eq!(r, [1.0, 2.0, 3.0]);
    }
}
