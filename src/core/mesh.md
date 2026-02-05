# Mesh.rs Mapping

| Rust Function | JS Function | JS Line Number | Description | Notes |
|---|---|---|---|---|
| `Mesh::new` | `Mesh.constructor` | 23 | Initializes a new empty Mesh. | |
| `Mesh::euler_characteristic` | `Mesh.eulerCharacteristic` | 38 | Computes the Euler characteristic (V - E + F). | |
| `Mesh::build` | `Mesh.build` | 51 | Constructs the mesh connectivity from a polygon soup. | |
| `Mesh::preallocate_elements` | `Mesh.preallocateElements` | 179 | Preallocates memory for mesh elements. | |
| `Mesh::has_isolated_vertices` | `Mesh.hasIsolatedVertices` | 224 | Checks for isolated vertices. | |
| `Mesh::has_isolated_faces` | `Mesh.hasIsolatedFaces` | 239 | Checks for isolated faces. | |
| `Mesh::has_non_manifold_vertices` | `Mesh.hasNonManifoldVertices` | 258 | Checks for non-manifold vertices. | |
| `Mesh::index_elements` | `Mesh.indexElements` | 283 | Assigns indices to all mesh elements. | |
| `Mesh::index_elements_global` | `indexElements` | 316 | Assigns unique indices to a list of elements. | Global JS function, but redundant with `index_elements` which populates indices internally |
| `Mesh::vertex_degree` | (inline logic) | | logic to compute vertex degree. | Rust helper for validation |
| `Mesh::is_isolated` | (inline logic) | | logic to check isolation. | Rust helper |
| `Mesh::on_boundary` | (inline logic) | | logic to check boundary status. | Rust helper |

## Iterators
| Rust Function | JS Function | Description | Notes |
|---|---|---|---|
| `vertex_adjacent_vertices` | `VertexVertexIterator` | Iterates over vertices adjacent to a vertex. | Implemented as Rust Iterator |
| `vertex_adjacent_edges` | `VertexEdgeIterator` | Iterates over edges adjacent to a vertex. | |
| `vertex_adjacent_faces` | `VertexFaceIterator` | Iterates over faces adjacent to a vertex. | |
| `vertex_adjacent_halfedges` | `VertexHalfedgeIterator` | Iterates over halfedges adjacent to a vertex. | |
| `vertex_adjacent_corners` | `VertexCornerIterator` | Iterates over corners adjacent to a vertex. | |
| `face_adjacent_vertices` | `FaceVertexIterator` | Iterates over vertices in a face. | |
| `face_adjacent_edges` | `FaceEdgeIterator` | Iterates over edges in a face. | |
| `face_adjacent_faces` | `FaceFaceIterator` | Iterates over faces adjacent to a face. | |
| `face_adjacent_halfedges` | `FaceHalfedgeIterator` | Iterates over halfedges in a face. | |
| `face_adjacent_corners` | `FaceCornerIterator` | Iterates over corners in a face. | |
| `boundary_adjacent_vertices` | (Implicit in JS) | Iterates over vertices in a boundary loop. | |
