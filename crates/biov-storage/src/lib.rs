//! Offline, portable copies of complete native packages. See `README.md` for
//! validation limits, content encoding and trusted-filesystem assumptions.
mod native;
mod store;
mod tree;
mod types;
pub use store::NativeStore;
pub use types::*;

#[cfg(test)]
mod tests;
