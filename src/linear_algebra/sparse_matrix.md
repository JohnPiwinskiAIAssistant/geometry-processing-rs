# SparseMatrix.rs Mapping

| Rust Function | JS Function | JS Line Number | Description | Notes |
|---|---|---|---|---|
| `SparseMatrix::from_triplet` | `SparseMatrix.fromTriplet` | 15 | Constructs sparse matrix from triplets. | |
| `SparseMatrix::get` | (indexing) | | Gets value at (row, col). | |
| `SparseMatrix::n_rows` | `SparseMatrix.nRows` | 33 | Returns number of rows. | |
| `SparseMatrix::n_cols` | `SparseMatrix.nCols` | 42 | Returns number of columns. | |
| `SparseMatrix::nnz` | `SparseMatrix.nnz` | 51 | Returns number of non-zero entries. | |
| `SparseMatrix::frobenius_norm` | `SparseMatrix.frobeniusNorm` | 60 | Computes Frobenius norm. | |
| `SparseMatrix::to_dense` | `SparseMatrix.toDense` | 74 | Converts to dense matrix. | |
| `SparseMatrix::plus` | `SparseMatrix.plus` | 84 | Adds two sparse matrices. | |
| `SparseMatrix::times_real` | `SparseMatrix.times` (scalar) | 99 | Multiplies by a real scalar. | |
| `SparseMatrix::times_dense` | `SparseMatrix.timesDense` | 114 | Multiplies by a dense matrix. | |
| `SparseMatrix::times_sparse` | `SparseMatrix.timesSparse` | 129 | Multiplies by another sparse matrix. | |
| `SparseMatrix::transpose` | `SparseMatrix.transpose` | 144 | Returns transpose. | |
| `SparseMatrix::invert_diagonal` | `SparseMatrix.invertDiagonal` | 159 | Inverts diagonal elements (for preconditioning). | |
| `SparseMatrix::identity` | `SparseMatrix.identity` | 175 | Creates identity sparse matrix. | |
| `SparseMatrix::diag` | `SparseMatrix.diag` | 186 | Creates diagonal sparse matrix from vector. | |
| `SparseMatrix::sub_matrix` | `SparseMatrix.subMatrix` | 197 | Extracts submatrix. | |
| `SparseMatrix::increment_by` | `SparseMatrix.incrementBy` | 211 | Adds another matrix in place. | |
| `SparseMatrix::scale_by` | `SparseMatrix.scaleBy` | 222 | Scales by scalar in place. | |
| `SparseMatrix::chol` | `SparseMatrix.chol` | 236 | Computes Cholesky decomposition. | Returns `Cholesky` solver object |
| `SparseMatrix::lu` | `SparseMatrix.lu` | 247 | Computes LU decomposition. | Returns `LU` solver object |
| `SparseMatrix::qr` | `SparseMatrix.qr` | 258 | Computes QR decomposition. | Returns `QR` solver object |

## Triplet.rs Mapping

| Rust Function | JS Function | JS Line Number | Description | Notes |
|---|---|---|---|---|
| `Triplet::new` | `Triplet` (constructor) | 10 | Initializes triplet list. | |
| `Triplet::add_entry` | `Triplet.addEntry` | 19 | Adds an entry (val, row, col). | |

## ComplexSparseMatrix.rs Mapping

| Rust Function | JS Function | JS Line Number | Description | Notes |
|---|---|---|---|---|
| `ComplexSparseMatrix::from_triplet` | `ComplexSparseMatrix.fromTriplet` | | Constructs from complex triplets. | |
| `ComplexSparseMatrix::n_rows` | `ComplexSparseMatrix.nRows` | | Returns number of rows. | |
| `ComplexSparseMatrix::n_cols` | `ComplexSparseMatrix.nCols` | | Returns number of columns. | |
| `ComplexSparseMatrix::nnz` | `ComplexSparseMatrix.nnz` | | Returns number of non-zero entries. | |
| `ComplexSparseMatrix::frobenius_norm` | `ComplexSparseMatrix.frobeniusNorm` | | Computes Frobenius norm. | |
| `ComplexSparseMatrix::plus` | `ComplexSparseMatrix.plus` | | Adds two complex sparse matrices. | |
| `ComplexSparseMatrix::minus` | `ComplexSparseMatrix.minus` | | Subtracts two complex sparse matrices. | |
| `ComplexSparseMatrix::scale_by` | `ComplexSparseMatrix.scaleBy` | | Scales by complex scalar in place. | |
| `ComplexSparseMatrix::to_dense` | `ComplexSparseMatrix.toDense` | | Converts to complex dense matrix. | |
| (Missing) | `ComplexSparseMatrix.chol` | | Cholesky decomposition. | **Missing in Rust port** |
| (Missing) | `ComplexSparseMatrix.lu` | | LU decomposition. | **Missing in Rust port** |
| (Missing) | `ComplexSparseMatrix.qr` | | QR decomposition. | **Missing in Rust port** |
