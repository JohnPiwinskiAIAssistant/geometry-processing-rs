# Poisson Problem

The Poisson problem solver computes scalar fields on surfaces by solving the Poisson equation $\Delta \phi = \rho$, where $\Delta$ is the Laplace-Beltrami operator and $\rho$ is a given source term.

```rust
pub struct ScalarPoissonProblem<'a> {
    pub geometry: &'a Geometry<'a>,
    pub a: SparseMatrix,      // Laplace matrix
    pub m: SparseMatrix,      // Mass matrix
    pub total_area: f64,
}
```

## Constructor

```rust
pub fn new(geometry: &'a Geometry<'a>) -> Self
```

Initializes the Poisson solver by precomputing:
- The Laplace matrix $L$ (cotan-Laplace operator)
- The mass matrix $M$ (diagonal matrix of vertex areas)
- The total surface area

These matrices are reused across multiple solves with different source terms.

## Solving the Poisson Equation

```rust
pub fn solve(&self, rho: &DenseMatrix) -> DenseMatrix
```

Solves the Poisson equation $L\phi = M(\bar{\rho} - \rho)$ for the scalar field $\phi$, where:
- $\rho$ is the input source term (one value per vertex)
- $\bar{\rho}$ is the mean value of $\rho$ weighted by vertex areas
- $L$ is the Laplace matrix
- $M$ is the mass matrix

### Algorithm

1. **Compute mean source**: $\bar{\rho} = \frac{\sum_i M_{ii} \rho_i}{\text{total area}}$

2. **Build right-hand side**: $b = M(\bar{\rho} - \rho)$
   - This ensures the solution exists (makes the system compatible)
   - The Laplacian has a null space (constant functions), so we need $\int \rho = 0$

3. **Solve linear system**: $L\phi = b$ using Cholesky factorization
   - The Laplace matrix is symmetric positive semi-definite
   - The solution is unique up to an additive constant

### Applications

The Poisson equation appears in many geometry processing tasks:
- **Surface reconstruction**: Given gradients or normals, reconstruct heights
- **Seamless texture synthesis**: Blend textures by solving for pixel values
- **Mesh smoothing**: Smooth a scalar field while preserving features
- **Heat diffusion**: Model heat flow on surfaces (related to geodesic distance)

### Mathematical Background

On a smooth surface, the Poisson equation is:
$$\Delta \phi = \rho$$

where $\Delta$ is the Laplace-Beltrami operator. In the discrete setting, we use the cotan-Laplace matrix as an approximation.

The compatibility condition $\int_S \rho \, dA = 0$ must hold for a solution to exist. This is why we subtract the mean $\bar{\rho}$ from $\rho$ before solving.
