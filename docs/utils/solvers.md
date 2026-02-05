# Solvers

The `Solvers` module provides numerical methods for solving eigenvalue problems and other computational tasks.

```rust
pub struct Solvers;
```

This is a utility struct with static methods for various numerical algorithms.

## Inverse Power Method

```rust
pub fn solve_inverse_power_method(a: &ComplexSparseMatrix) -> ComplexDenseMatrix
```

Computes the eigenvector corresponding to the **smallest eigenvalue** of a complex sparse matrix using the inverse power method.

### Algorithm

The inverse power method is an iterative algorithm:

1. **Initialize**: Start with a random vector $x_0$
2. **Iterate**: For $k = 1, 2, \ldots, 200$:
   - Solve $Ay_k = x_{k-1}$ for $y_k$ (using sparse LU factorization)
   - Normalize: $x_k = y_k / \|y_k\|$
3. **Return**: The final $x_k$ approximates the eigenvector

### Why It Works

If $A$ has eigenvalues $\lambda_1 \leq \lambda_2 \leq \cdots \leq \lambda_n$ with corresponding eigenvectors $v_1, v_2, \ldots, v_n$, then:
- Multiplying by $A^{-1}$ amplifies components in the direction of $v_1$ (smallest eigenvalue)
- After many iterations, $x_k \approx v_1$

This is used in spectral conformal parameterization to find the conformal map that minimizes distortion.

### Implementation Details

- Uses `faer`'s sparse LU factorization (`sp_lu`) for efficiency
- Runs for 200 iterations (empirically sufficient for convergence)
- Normalizes at each step to prevent numerical overflow/underflow

## Residual Computation

```rust
pub fn residual(a: &ComplexSparseMatrix, x: &ComplexDenseMatrix) -> f64
```

Computes the residual $\|Ax\|_2$ for a given matrix $A$ and vector $x$. This is useful for:
- Checking convergence of iterative methods
- Verifying that $x$ is close to a null vector (when $x$ is an approximate eigenvector for eigenvalue $\approx 0$)

The implementation manually performs sparse matrix-vector multiplication to compute $Ax$, then returns its 2-norm.

## Matrix Inversion

```rust
pub fn invert_2x2(m: &DenseMatrix) -> DenseMatrix
```

Inverts a $2 \times 2$ matrix using the explicit formula:
$$A^{-1} = \frac{1}{\det(A)} \begin{bmatrix} d & -b \\ -c & a \end{bmatrix}$$

where $A = \begin{bmatrix} a & b \\ c & d \end{bmatrix}$ and $\det(A) = ad - bc$.

This is used in various geometry processing algorithms that require local 2D transformations.
