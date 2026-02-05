# Dense Matrix

The `DenseMatrix` and `ComplexDenseMatrix` structs provide dense matrix storage and operations. These are wrappers around the `faer` crate's matrix types.

## DenseMatrix

```rust
pub struct DenseMatrix {
    pub inner: Mat<f64>,
}
```

A dense matrix of real numbers stored in column-major format.

### Constructors

```rust
pub fn zeros(m: usize, n: usize) -> Self
```
Creates an $m \times n$ matrix filled with zeros.

---

```rust
pub fn identity(m: usize, n: usize) -> Self
```
Creates an $m \times n$ identity matrix (1s on diagonal, 0s elsewhere).

---

```rust
pub fn ones(m: usize, n: usize) -> Self
```
Creates an $m \times n$ matrix filled with ones.

---

```rust
pub fn constant(x: f64, m: usize, n: usize) -> Self
```
Creates an $m \times n$ matrix filled with the constant value `x`.

---

```rust
pub fn random(m: usize, n: usize) -> Self
```
Creates an $m \times n$ matrix with random entries in $[0, 1)$.

### Properties

```rust
pub fn n_rows(&self) -> usize
pub fn n_cols(&self) -> usize
```
Returns the number of rows and columns.

---

```rust
pub fn transpose(&self) -> Self
```
Returns the transpose $A^T$.

---

```rust
pub fn rank(&self) -> usize
```
Computes the rank of the matrix (number of linearly independent rows/columns).

### Norms

```rust
pub fn norm(&self, n: usize) -> f64
```
Computes matrix norms:
- `n = 0`: Infinity norm (maximum absolute row sum)
- `n = 1`: 1-norm (maximum absolute column sum)  
- `n = 2` (or other): Frobenius norm $\sqrt{\sum_{i,j} |a_{ij}|^2}$

### Operations

```rust
pub fn get(&self, i: usize, j: usize) -> f64
pub fn set(&mut self, x: f64, i: usize, j: usize)
```
Access and modify individual elements.

---

```rust
pub fn plus(&self, b: &DenseMatrix) -> Self      // A + B
pub fn minus(&self, b: &DenseMatrix) -> Self     // A - B
pub fn times_dense(&self, b: &DenseMatrix) -> Self  // A * B
pub fn scale_by(&mut self, s: f64)               // A *= s
```

Arithmetic operations on matrices.

---

```rust
pub fn increment_by(&mut self, b: &DenseMatrix)  // A += B
pub fn decrement_by(&mut self, b: &DenseMatrix)  // A -= B
```

In-place operations.

---

```rust
pub fn negate(&mut self)
pub fn negated(&self) -> Self
```

Negation (in-place and returning new matrix).

## ComplexDenseMatrix

```rust
pub struct ComplexDenseMatrix {
    pub mat: Vec<Complex64>,
    pub rows: usize,
    pub cols: usize,
}
```

A dense matrix of complex numbers. Unlike `DenseMatrix`, this uses a custom storage format with a `Vec<Complex64>`.

### Constructors

```rust
pub fn zeros(m: usize, n: usize) -> Self
pub fn identity(m: usize, n: usize) -> Self
pub fn ones(m: usize, n: usize) -> Self
pub fn constant(x: Complex, m: usize, n: usize) -> Self
pub fn random(m: usize, n: usize) -> Self
```

Similar to `DenseMatrix` but for complex values.

### Operations

The API mirrors `DenseMatrix` but operates on `Complex` values:

```rust
pub fn get(&self, i: usize, j: usize) -> Complex
pub fn set(&mut self, x: Complex, i: usize, j: usize)
pub fn conjugate(&self) -> Self
pub fn transpose(&self) -> Self
pub fn norm(&self, n: usize) -> f64
```

The `conjugate()` method returns the complex conjugate of each element, and `transpose()` returns the Hermitian transpose (conjugate transpose).
