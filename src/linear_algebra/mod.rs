pub mod sparse_matrix;
pub mod traits;

pub type Vector = faer::Mat<f64>;
pub type PolygonSoup = crate::core::mesh::PolygonSoup;
pub type Complex = num_complex::Complex64;
pub type DenseMatrix = faer::Mat<f64>;
pub type ComplexDenseMatrix = faer::Mat<Complex>;
pub type SparseMatrix<S> = faer::sparse::SparseColMat<usize, S>;
pub use sparse_matrix::{Triplet, Cholesky, LU, QR};
