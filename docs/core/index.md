# Core

This module contains the halfedge mesh data structure and geometric operations. It forms the foundation of the geometry processing framework.

## Mesh Data Structure

The halfedge mesh is a powerful representation that makes it easy to traverse and query mesh connectivity.

### [Mesh](./mesh.md)
The main mesh class containing vertices, edges, faces, halfedges, corners, and boundaries.

### [Vertex](./vertex.md)
Represents a vertex in the mesh.

### [Edge](./edge.md)
Represents an undirected edge connecting two vertices.

### [Face](./face.md)
Represents a face (typically a triangle) in the mesh.

### [Halfedge](./halfedge.md)
Represents a directed edge with connectivity information.

### [Corner](./corner.md)
Represents a corner (vertex of a face) with associated angle information.

## Geometric Operations

### [Geometry](./geometry.md)
Provides geometric measurements and operators on the mesh (edge lengths, face areas, curvatures, Laplacian matrices, etc.).

### [DEC](./dec.md)
Discrete Exterior Calculus operators (Hodge star, exterior derivative, etc.).
