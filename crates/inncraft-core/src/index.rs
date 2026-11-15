//! Index types.
//!
//! Local indices address cells, faces and matrix entries inside one rank or
//! GPU and are 32 bit, which halves the index traffic of the kernels against
//! 64 bit. At the design target of 100+ million cells (D6) a single partition
//! stays far below 2^31 entries. Global indices, used for I/O, restart with a
//! different decomposition (D42) and particle identities (D47), are 64 bit.

/// Index of a cell, face or matrix entry inside one partition.
pub type LocalIndex = u32;

/// Index of a cell or face in the whole, undecomposed mesh.
pub type GlobalIndex = u64;

/// Largest local index in use. The top bit of a [`LocalIndex`] is reserved
/// for flags packed into connectivity arrays (for example the orientation of
/// a face seen from a cell).
pub const MAX_LOCAL_INDEX: LocalIndex = (1 << 31) - 1;

/// Flag bit of a packed [`LocalIndex`].
pub const FLAG_BIT: LocalIndex = 1 << 31;
