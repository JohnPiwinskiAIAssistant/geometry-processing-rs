# Geometry

The `Geometry` struct provides geometric measurements and operators on a mesh. It stores vertex positions and provides methods to compute edge lengths, face areas, curvatures, and discrete differential operators.

```rust
pub struct Geometry<'a> {
    pub mesh: &'a Mesh,
    pub positions: Vec<Vector>,
}
```

## Construction

```rust
pub fn new(mesh: &'a Mesh, positions: Vec<Vector>, normalize_positions: bool) -> Self
```

Creates a new geometry from a mesh and vertex positions. If `normalize_positions` is true, the mesh is centered at the origin and scaled to unit radius.

---

```rust
pub fn normalize(&mut self, rescale: bool)
```

Centers the mesh at the origin by translating all positions by the negative center of mass. If `rescale` is true, also scales positions so the furthest vertex is at distance 1.

## Basic Measurements

```rust
pub fn vector(&self, h_idx: usize) -> Vector
```

Returns the vector along a halfedge from its tail vertex to its head vertex.

---

```rust
pub fn length(&self, e_idx: usize) -> f64
pub fn midpoint(&self, e_idx: usize) -> Vector
pub fn mean_edge_length(&self) -> f64
```

Edge measurements: length, midpoint, and average length across all edges.

---

```rust
pub fn area(&self, f: &Face) -> f64
pub fn total_area(&self) -> f64
```

Face area (using cross product formula) and total surface area.

---

```rust
pub fn face_normal(&self, f: &Face) -> Option<Vector>
pub fn centroid(&self, f: &Face) -> Vector
pub fn circumcenter(&self, f: &Face) -> Vector
pub fn orthonormal_bases(&self, f: &Face) -> [Vector; 2]
```

Face properties: normal vector, centroid (average of vertices), circumcenter (center of circumscribed circle), and orthonormal basis vectors tangent to the face.

## Angles and Cotangents

```rust
pub fn angle(&self, c: &Corner) -> f64
```

Computes the interior angle at a corner (in radians), clamped to $[0, \pi]$.

---

```rust
pub fn cotan(&self, h_idx: usize) -> f64
```

Computes the cotangent of the angle opposite to a halfedge. This is a fundamental quantity in discrete differential geometry, appearing in the Laplace-Beltrami operator. Returns 0 for boundary halfedges.

---

```rust
pub fn dihedral_angle(&self, h_idx: usize) -> f64
```

Computes the signed dihedral angle between two adjacent faces sharing a halfedge.

## Dual Areas

```rust
pub fn barycentric_dual_area(&self, v: &Vertex) -> f64
```

Computes the barycentric dual area: sum of 1/3 of the area of each adjacent face. This is used in the mass matrix.

---

```rust
pub fn circumcentric_dual_area(&self, v: &Vertex) -> f64
```

Computes the circumcentric dual area using the formula involving cotangents. This is more accurate for well-shaped triangles.

## Vertex Normals

The geometry provides multiple methods for computing vertex normals, each with different weighting schemes:

```rust
pub fn vertex_normal_equally_weighted(&self, v: &Vertex) -> Vector
```
Average of adjacent face normals (uniform weights).

---

```rust
pub fn vertex_normal_area_weighted(&self, v: &Vertex) -> Vector
```
Weighted by face areas.

---

```rust
pub fn vertex_normal_angle_weighted(&self, v: &Vertex) -> Vector
```
Weighted by corner angles (tip angles).

---

```rust
pub fn vertex_normal_gauss_curvature(&self, v: &Vertex) -> Vector
```
Based on Gauss curvature formula using dihedral angles.

---

```rust
pub fn vertex_normal_mean_curvature(&self, v: &Vertex) -> Vector
```
Based on mean curvature formula (same as area gradient method).

---

```rust
pub fn vertex_normal_sphere_inscribed(&self, v: &Vertex) -> Vector
```
Based on inscribed sphere method.

## Curvatures

```rust
pub fn angle_defect(&self, v: &Vertex) -> f64
```

Computes the angle defect: $2\pi - \sum \theta_i$ for interior vertices, $\pi - \sum \theta_i$ for boundary vertices. This equals the integrated Gaussian curvature.

---

```rust
pub fn scalar_gauss_curvature(&self, v: &Vertex) -> f64
pub fn total_angle_defect(&self) -> f64
```

Gaussian curvature (equals angle defect) and total angle defect (equals $2\pi \chi$ by Gauss-Bonnet).

---

```rust
pub fn scalar_mean_curvature(&self, v: &Vertex) -> f64
```

Integrated mean curvature using dihedral angles and edge lengths.

---

```rust
pub fn principal_curvatures(&self, v: &Vertex) -> [f64; 2]
```

Computes the minimum and maximum principal curvatures $\kappa_1, \kappa_2$ from the mean and Gaussian curvatures.

## Discrete Operators

```rust
pub fn laplace_matrix(&self) -> SparseMatrix
```

Builds the cotan-Laplace matrix $L$, where $L_{ij} = -\frac{1}{2}(\cot \alpha + \cot \beta)$ for edge $(i,j)$ and $L_{ii} = -\sum_{j \neq i} L_{ij}$. This is a discrete approximation of the Laplace-Beltrami operator.

---

```rust
pub fn mass_matrix(&self) -> SparseMatrix
```

Builds the diagonal mass matrix $M$ with $M_{ii}$ equal to the barycentric dual area of vertex $i$.

---

```rust
pub fn complex_laplace_matrix(&self) -> ComplexSparseMatrix
```

Complex version of the Laplace matrix, used in algorithms that work with complex-valued functions (e.g., conformal parameterization).
