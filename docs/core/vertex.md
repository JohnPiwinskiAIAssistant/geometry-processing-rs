# Vertex.rs Mapping

| Rust Function | JS Function | JS Line Number | Description | Notes |
|---|---|---|---|---|
| `Vertex::new` | `Vertex.constructor` | 10 | Initializes a new Vertex instance. | |
| `Vertex::degree` | `Vertex.degree` | 19 | Counts the number of edges adjacent to this vertex. | Takes `&Mesh` |
| `Vertex::is_isolated` | `Vertex.isIsolated` | 30 | Checks whether this vertex is isolated (has no incident halfedge). | |
| `Vertex::on_boundary` | `Vertex.onBoundary` | 39 | Checks whether this vertex lies on a boundary. | Takes `&Mesh` |
| `Vertex::adjacent_vertices` | `Vertex.adjacentVertices` | 55 | Iterates over vertices neighboring this vertex. | Takes `&Mesh` |
| `Vertex::adjacent_edges` | `Vertex.adjacentEdges` | 69 | Iterates over edges adjacent to this vertex. | Takes `&Mesh` |
| `Vertex::adjacent_faces` | `Vertex.adjacentFaces` | 83 | Iterates over faces adjacent to this vertex. | Takes `&Mesh` |
| `Vertex::adjacent_halfedges` | `Vertex.adjacentHalfedges` | 97 | Iterates over outgoing halfedges from this vertex. | Takes `&Mesh` |
| `Vertex::adjacent_corners` | `Vertex.adjacentCorners` | 111 | Iterates over corners adjacent to this vertex. | Takes `&Mesh` |
| (Not implemented) | `Vertex.toString` | 122 | Returns a string representation (index). | Implemented via `Debug` trait |
