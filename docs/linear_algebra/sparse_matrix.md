# Sparse Matrix

The `SparseMatrix` and `ComplexSparseMatrix` structs provide sparse matrix storage in Compressed Column Storage (CCS) format, along with factorizations for solving linear systems.

## SparseMatrix

```rust
pub struct SparseMatrix {
    pub mat: Arc<SparseColMat<usize, f64>>,
}
```

A sparse matrix of real numbers. Uses `faer`'s sparse matrix implementation internally.

### Construction

```rust
pub fn from_triplet(triplet: Triplet) -> Self
```

Constructs a sparse matrix from a `Triplet` (list of $(i, j, value)$ entries). This is the primary way to build sparse matrices.

### Triplet

```rust
pub struct Triplet {
    pub m: usize,
    pub n: usize,
    pub entries: Vec<(usize, usize, f64)>,
}
```

A triplet stores matrix entries as a list of $(row, col, value)$ tuples before converting to CCS format.

```rust
impl Triplet {
    pub fn new(m: usize, n: usize) -> Self
    pub fn add_entry(&mut self, x: f64, i: usize, j: usize)
}
```

### Properties

```rust
pub fn n_rows(&self) -> usize
pub fn n_cols(&self) -> usize
pub fn nnz(&self) -> usize  // Number of non-zero entries
```

---

```rust
pub fn get(&self, i: usize, j: usize) -> f64
```

Retrieves an element. Returns 0.0 if the entry is not stored.

### Operations

```rust
pub fn transpose(&self) -> Self
pub fn plus(&self, b: &SparseMatrix) -> Self
pub fn minus(&self, b: &SparseMatrix) -> Self
pub fn times_sparse(&self, b: &SparseMatrix) -> Self
pub fn times_dense(&self, b: &DenseMatrix) -> DenseMatrix
pub fn scale_by(&mut self, s: f64)
```

Arithmetic operations. Note that `times_sparse` performs sparse-sparse multiplication, while `times_dense` performs sparse-dense multiplication.

---

```rust
pub fn to_dense(&self) -> DenseMatrix
```

Converts the sparse matrix to a dense matrix. Useful for debugging or when dense operations are needed.

### Linear System Solvers

The sparse matrix provides factorizations for solving linear systems $Ax = b$:

```rust
pub struct Cholesky { /* ... */ }
pub struct LU { /* ... */ }
pub struct QR { /* ... */ }
```

Each factorization is constructed from a sparse matrix and can then solve systems efficiently:

```rust
impl SparseMatrix {
    pub fn chol(&self) -> Cholesky
    pub fn lu(&self) -> LU
    pub fn qr(&self) -> QR
}

impl Cholesky {
    pub fn solve(&self, b: &DenseMatrix) -> DenseMatrix
}
// Similar for LU and QR
```

**Cholesky** factorization is used for symmetric positive-definite matrices (common in geometry processing). **LU** handles general square matrices. **QR** is useful for least-squares problems.

## ComplexSparseMatrix

```rust
pub struct ComplexSparseMatrix {
    pub mat: Arc<SparseColMat<usize, Complex64>>,
}
```

A sparse matrix of complex numbers, with an API mirroring `SparseMatrix`.

### ComplexTriplet

```rust
pub struct ComplexTriplet {
    pub m: usize,
    pub n: usize,
    pub entries: Vec<(usize, usize, Complex)>,
}
```

Used to construct `ComplexSparseMatrix` from a list of entries.

```rust
impl ComplexTriplet {
    pub fn new(m: usize, n: usize) -> Self
    pub fn add_entry(&mut self, x: Complex, i: usize, j: usize)
}
```

### Operations

```rust
pub fn conjugate(&self) -> Self
pub fn transpose(&self) -> Self
pub fn hermitian(&self) -> Self  // Conjugate transpose
pub fn frobenius_norm(&self) -> f64
```

The `hermitian()` method returns $A^H = \overline{A^T}$, the conjugate transpose, which is commonly used in complex linear algebra.

### Complex Linear Solvers

```rust
pub struct ComplexLU { /* ... */ }

impl ComplexSparseMatrix {
    pub fn lu(&self) -> ComplexLU
}

impl ComplexLU {
    pub fn solve(&self, b: &ComplexDenseMatrix) -> ComplexDenseMatrix
}
```

Currently, only LU factorization is implemented for complex sparse matrices.
