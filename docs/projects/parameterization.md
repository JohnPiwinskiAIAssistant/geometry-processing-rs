# Parameterization

This module implements algorithms for flattening 3D surfaces to 2D while preserving certain geometric properties. Parameterization is fundamental for texture mapping, remeshing, and other geometry processing tasks.

## Spectral Conformal Parameterization

```rust
pub struct SpectralConformalParameterization<'a> {
    pub geometry: &'a Geometry<'a>,
}
```

This algorithm computes a conformal (angle-preserving) parameterization by finding the eigenvector corresponding to the smallest eigenvalue of the conformal energy matrix.

### Constructor

```rust
pub fn new(geometry: &'a Geometry<'a>) -> Self
```

Creates a new spectral conformal parameterization solver for the given geometry.

### Conformal Energy Matrix

```rust
pub fn build_conformal_energy(&self) -> ComplexSparseMatrix
```

Builds the complex conformal energy matrix $E_C = E_D - A$, where:
- $E_D$ is the Dirichlet energy (complex Laplacian scaled by 0.5)
- $A$ is the area term encoding boundary constraints

The matrix is constructed as follows:
1. Start with $E_D = \frac{1}{2}L$, where $L$ is the complex Laplace matrix
2. For each boundary halfedge $(i, j)$, add:
   - $A_{ij} = \frac{i}{4}$ (where $i = \sqrt{-1}$)
   - $A_{ji} = -\frac{i}{4}$
3. Return $E_C = E_D - A$

This formulation comes from the discrete conformal parameterization theory, where the conformal energy measures deviation from conformality.

### Flattening

```rust
pub fn flatten(&self) -> Vec<Vector>
```

Computes the 2D parameterization by:
1. Building the conformal energy matrix $E_C$
2. Finding the eigenvector $z$ corresponding to the smallest eigenvalue using inverse power iteration
3. Interpreting $z$ as complex coordinates: $u_i = \text{Re}(z_i)$, $v_i = \text{Im}(z_i)$
4. Normalizing the result (centering and scaling to unit radius)

The resulting parameterization minimizes conformal distortion while satisfying boundary constraints.

### Normalization

```rust
fn normalize_flattening(&self, flattening: &mut [Vector])
```

Centers the flattening at the origin and scales it to unit radius. This is done for numerical stability and to provide a canonical representation.

## Boundary First Flattening

```rust
pub struct BoundaryFirstFlattening<'a> {
    pub geometry: &'a Geometry<'a>,
    pub n_v: usize,        // Number of vertices
    pub n_i: usize,        // Number of interior vertices
    pub n_b: usize,        // Number of boundary vertices
    pub vertex_index: HashMap<usize, usize>,
    pub b_vertex_index: HashMap<usize, usize>,
    pub k_gaussian: DenseMatrix,   // Gaussian curvatures
    pub k_geodesic: DenseMatrix,   // Geodesic curvatures
    pub l_lengths: DenseMatrix,    // Edge lengths
    pub a_ii: SparseMatrix,        // Interior-interior block
    pub a_ib: SparseMatrix,        // Interior-boundary block
    pub a_bb: SparseMatrix,        // Boundary-boundary block
    pub a_full: SparseMatrix,      // Full system matrix
}
```

Boundary First Flattening (BFF) is a more sophisticated algorithm that:
1. First solves for optimal boundary positions
2. Then solves for interior vertex positions

This two-stage approach allows for better control over distortion.

### Constructor

```rust
pub fn new(geometry: &'a Geometry<'a>) -> Self
```

Initializes the BFF solver by:
- Separating interior and boundary vertices
- Computing Gaussian curvatures at vertices
- Computing geodesic curvatures along boundary edges
- Building the system matrices

### Flattening

```rust
pub fn flatten(&self, target_boundary: Option<Vec<f64>>) -> Vec<Vector>
```

Computes the parameterization in two stages:

**Stage 1: Boundary Parameterization**
- If `target_boundary` is provided, use those target angles
- Otherwise, solve for optimal boundary angles that minimize distortion
- Map boundary vertices to a circle or other target shape

**Stage 2: Interior Parameterization**
- With boundary positions fixed, solve a Laplace equation for interior positions
- This ensures the interior mapping is harmonic (minimizes Dirichlet energy)

The result is a low-distortion parameterization with explicit control over the boundary shape.

### Implementation Details

The BFF algorithm uses:
- **Gaussian curvature** $K$ at interior vertices (angle defect)
- **Geodesic curvature** $\kappa_g$ along boundary edges
- **Cotan-Laplace operator** for the interior solve

The boundary solve involves finding scale factors $u_i$ such that the flattened boundary has the desired total curvature, which is formulated as a sparse linear system.
