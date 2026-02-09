pub mod sparse_matrix;

pub type Vector = faer::Mat<f64>;
pub type PolygonSoup = crate::core::mesh::PolygonSoup;
pub type Complex = num_complex::Complex64;
pub type DenseMatrix = faer::Mat<f64>;
pub type ComplexDenseMatrix = faer::Mat<Complex>;
pub type SparseMatrix = faer::sparse::SparseColMat<usize, f64>;
pub type ComplexSparseMatrix = faer::sparse::SparseColMat<usize, Complex>;
pub use sparse_matrix::{Triplet, ComplexTriplet, Cholesky, LU, QR};
