# DenseMatrix.rs Mapping

| Rust Function | JS Function | JS Line Number | Description | Notes |
|---|---|---|---|---|
| `DenseMatrix::from_inner` | (No match) | | Wraps a `faer::Mat`. | Rust internal |
| `DenseMatrix::zeros` | `DenseMatrix.zeros` | 15 | Creates an m x n matrix of zeros. | |
| `DenseMatrix::identity` | `DenseMatrix.identity` | 24 | Creates an m x n identity matrix. | |
| `DenseMatrix::ones` | `DenseMatrix.ones` | 33 | Creates an m x n matrix of ones. | |
| `DenseMatrix::constant` | `DenseMatrix.constant` | 42 | Creates an m x n matrix with constant value. | |
| `DenseMatrix::random` | `DenseMatrix.random` | 51 | Creates an m x n matrix with random values. | |
| `DenseMatrix::transpose` | `DenseMatrix.transpose` | 60 | Returns the transpose of the matrix. | |
| `DenseMatrix::n_rows` | `DenseMatrix.nRows` | 69 | Returns the number of rows. | |
| `DenseMatrix::n_cols` | `DenseMatrix.nCols` | 78 | Returns the number of columns. | |
| `DenseMatrix::norm` | `DenseMatrix.norm` | 87 | Computes vector norms (1, infinity, or Frobenius). | |
| `DenseMatrix::rank` | `DenseMatrix.rank` | 100 | Computes the rank of the matrix. | Currently a stub in Rust? |
| `DenseMatrix::sum` | `DenseMatrix.sum` | 110 | Computes the sum of all elements. | |
| `DenseMatrix::sub_matrix` | `DenseMatrix.subMatrix` | 120 | Extracts a submatrix. | |
| `DenseMatrix::increment_by` | `DenseMatrix.incrementBy` | 131 | Adds another matrix in place. | |
| `DenseMatrix::decrement_by` | `DenseMatrix.decrementBy` | 142 | Subtracts another matrix in place. | |
| `DenseMatrix::scale_by` | `DenseMatrix.scaleBy` | 153 | Scales the matrix by a scalar in place. | |
| `DenseMatrix::plus` | `DenseMatrix.plus` | 164 | Returns the sum of two matrices. | |
| `DenseMatrix::minus` | `DenseMatrix.minus` | 175 | Returns the difference of two matrices. | |
| `DenseMatrix::times_real` | `DenseMatrix.times` (scalar) | 186 | Multiplies by a real scalar. | |
| `DenseMatrix::times_dense` | `DenseMatrix.times` (matrix) | 186 | Multiplies by another matrix. | |
| `DenseMatrix::negated` | `DenseMatrix.negated` | 197 | Returns the negation of the matrix. | |
| `DenseMatrix::get` | (indexing) | | Gets value at (i, j). | |
| `DenseMatrix::set` | (indexing) | | Sets value at (i, j). | |
| `DenseMatrix::hcat` | `DenseMatrix.hcat` | 210 | Concatenates matrices horizontally. | |
| `DenseMatrix::vcat` | `DenseMatrix.vcat` | 225 | Concatenates matrices vertically. | |

## ComplexDenseMatrix.rs Mapping

| Rust Function | JS Function | JS Line Number | Description | Notes |
|---|---|---|---|---|
| `ComplexDenseMatrix::zeros` | `ComplexDenseMatrix.zeros` | | Creates an m x n complex matrix of zeros. | Custom impl |
| `ComplexDenseMatrix::ones` | `ComplexDenseMatrix.ones` | | Creates an m x n complex matrix of ones. | |
| `ComplexDenseMatrix::random` | `ComplexDenseMatrix.random` | | Creates an m x n complex matrix of random values. | |
| `ComplexDenseMatrix::n_rows` | `ComplexDenseMatrix.nRows` | | Returns number of rows. | |
| `ComplexDenseMatrix::n_cols` | `ComplexDenseMatrix.nCols` | | Returns number of columns. | |
| `ComplexDenseMatrix::transpose` | `ComplexDenseMatrix.transpose` | | Returns transpose. | |
| `ComplexDenseMatrix::conjugate` | `ComplexDenseMatrix.conjugate` | | Returns conjugate. | |
| `ComplexDenseMatrix::get` | (indexing) | | Gets complex value at (i, j). | |
| `ComplexDenseMatrix::set` | (indexing) | | Sets complex value at (i, j). | |
| `ComplexDenseMatrix::sum` | `ComplexDenseMatrix.sum` | | Computes sum of elements. | |
| `ComplexDenseMatrix::norm` | `ComplexDenseMatrix.norm` | | Computes Frobenius norm. | Only norm=2 supported currently? |
| `ComplexDenseMatrix::decrement_by` | `ComplexDenseMatrix.decrementBy` | | Subtracts in place. | |
| `ComplexDenseMatrix::scale_by` | `ComplexDenseMatrix.scaleBy` | | Scales in place. | |
| `ComplexDenseMatrix::times_dense` | `ComplexDenseMatrix.times` | | Matrix multiplication. | O(N^3) implementation |
| `ComplexDenseMatrix::minus` | `ComplexDenseMatrix.minus` | | Subtracts matrices. | |
| `ComplexDenseMatrix::times_complex` | `ComplexDenseMatrix.times` (scalar) | | Multiplies by complex scalar. | |
