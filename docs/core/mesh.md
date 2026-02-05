# Mesh

The `Mesh` struct represents a halfedge mesh data structure, which provides efficient access to mesh connectivity and topology.

## Mesh Structure

```rust
pub struct Mesh {
    pub vertices: Vec<Vertex>,
    pub edges: Vec<Edge>,
    pub faces: Vec<Face>,
    pub corners: Vec<Corner>,
    pub halfedges: Vec<Halfedge>,
    pub boundaries: Vec<Face>,
    pub generators: Vec<Vec<usize>>,
}
```

The mesh stores all elements in vectors indexed by their unique IDs. This allows for O(1) access to any element.

- **vertices**: All vertices in the mesh
- **edges**: All undirected edges
- **faces**: All interior faces (typically triangles)
- **corners**: All corners (one per vertex per face)
- **halfedges**: All directed halfedges (two per edge)
- **boundaries**: Boundary loops represented as special faces
- **generators**: Homology generators for non-trivial topology

## PolygonSoup

```rust
pub struct PolygonSoup {
    pub v: Vec<Vector>,
    pub f: Vec<usize>,
}
```

A simple representation of mesh geometry before building the halfedge structure:
- `v`: Vertex positions
- `f`: Face indices (flattened, 3 indices per triangle)

## Construction

```rust
pub fn new() -> Self
```

Creates an empty mesh.

---

```rust
pub fn build(&mut self, polygon_soup: &PolygonSoup) -> bool
```

Builds the halfedge mesh from a `PolygonSoup`. This method:
1. Creates vertices, edges, faces, and halfedges
2. Establishes connectivity (next, prev, twin pointers)
3. Identifies boundary loops
4. Detects non-manifold edges (returns `false` if found)

The algorithm ensures that each edge has at most 2 incident faces (manifold condition).

## Topology

```rust
pub fn euler_characteristic(&self) -> i32
```

Computes the Euler characteristic $\chi = V - E + F$, where $V$ is the number of vertices, $E$ is the number of edges, and $F$ is the number of faces. For a closed surface:
- Sphere: $\chi = 2$
- Torus: $\chi = 0$
- Surface with $g$ handles: $\chi = 2 - 2g$

## Connectivity Queries

The mesh provides various methods to traverse connectivity:

```rust
pub fn vertex_adjacent_vertices(&self, v_idx: usize) -> Vec<usize>
pub fn vertex_adjacent_edges(&self, v_idx: usize) -> Vec<usize>
pub fn vertex_adjacent_faces(&self, v_idx: usize) -> Vec<usize>
pub fn vertex_adjacent_halfedges(&self, v_idx: usize) -> Vec<usize>
pub fn vertex_adjacent_corners(&self, v_idx: usize) -> Vec<usize>
```

These methods return the indices of elements adjacent to a given vertex.

---

```rust
pub fn face_adjacent_vertices(&self, f_idx: usize) -> Vec<usize>
pub fn face_adjacent_edges(&self, f_idx: usize) -> Vec<usize>
pub fn face_adjacent_faces(&self, f_idx: usize) -> Vec<usize>
pub fn face_adjacent_halfedges(&self, f_idx: usize, ccw: bool) -> Vec<usize>
pub fn face_adjacent_corners(&self, f_idx: usize) -> Vec<usize>
```

These methods return elements adjacent to a given face. The `ccw` parameter controls the order of halfedges (counter-clockwise if `true`).

---

```rust
pub fn edge_adjacent_vertices(&self, e_idx: usize) -> Vec<usize>
pub fn edge_adjacent_faces(&self, e_idx: usize) -> Vec<usize>
pub fn edge_adjacent_halfedges(&self, e_idx: usize) -> Vec<usize>
```

These methods return elements adjacent to a given edge.

## Boundary Detection

```rust
pub fn is_boundary_vertex(&self, v_idx: usize) -> bool
pub fn is_boundary_edge(&self, e_idx: usize) -> bool
```

Checks if a vertex or edge lies on the mesh boundary.

## Implementation Notes

The halfedge data structure is particularly efficient for:
- **Local traversal**: Walking around a vertex or face
- **Adjacency queries**: Finding neighbors of any element
- **Mesh editing**: Splitting edges, collapsing vertices, etc.

Each halfedge stores pointers to its `next`, `prev`, `twin`, associated `vertex`, `edge`, and `face`, enabling efficient traversal in any direction.
