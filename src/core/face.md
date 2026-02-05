# Face.rs Mapping

| Rust Function | JS Function | JS Line Number | Description | Notes |
|---|---|---|---|---|
| `Face::new` | `Face.constructor` | 9 | Initializes a new Face instance. | |
| `Face::is_boundary_loop` | `Face.isBoundaryLoop` | 21 | Checks whether this face is a boundary loop. | Takes `&Mesh` |
| `Face::adjacent_vertices` | `Face.adjacentVertices` | 36 | Iterates over the vertices in this face. | Takes `&Mesh` |
| `Face::adjacent_edges` | `Face.adjacentEdges` | 50 | Iterates over the edges in this face. | Takes `&Mesh` |
| `Face::adjacent_faces` | `Face.adjacentFaces` | 64 | Iterates over the faces neighboring this face. | Takes `&Mesh` |
| `Face::adjacent_halfedges` | `Face.adjacentHalfedges` | 78 | Iterates over the halfedges in this face. | Takes `&Mesh` |
| `Face::adjacent_corners` | `Face.adjacentCorners` | 92 | Iterates over the corners in this face. | Takes `&Mesh` |
| (Not implemented) | `Face.toString` | 103 | Returns a string representation (index). | Implemented via `Debug` trait |
