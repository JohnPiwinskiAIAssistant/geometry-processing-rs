# Geometry Processing Rust Port - Documentation

This directory contains comprehensive documentation for the Rust port of the `geometry-processing-js` framework.

## Quick Start

- **[Main Index](./index.md)**: Overview of all modules
- **[Linear Algebra](./linear_algebra/index.md)**: Vectors, matrices, and numerical types
- **[Core](./core/index.md)**: Halfedge mesh data structure and geometric operations
- **[Projects](./projects/index.md)**: Geometry processing algorithms
- **[Utils](./utils/index.md)**: Utility functions and helpers

## Documentation Structure

Each module has:
- An **index.md** file providing an overview
- Individual **.md** files for each source file, with:
  - Rust code blocks showing struct/trait/function signatures
  - Detailed explanations of functionality
  - Mathematical background where relevant
  - Usage examples and implementation notes

## Key Modules

### Linear Algebra
Foundation for numerical computation:
- `Complex`: Complex numbers
- `Vector`: 3D vectors
- `DenseMatrix`, `SparseMatrix`: Matrix types with solvers

### Core
Halfedge mesh representation:
- `Mesh`: Connectivity and topology
- `Geometry`: Measurements and discrete operators
- `Vertex`, `Edge`, `Face`, `Halfedge`, `Corner`: Mesh elements

### Projects
Geometry processing algorithms:
- `Parameterization`: Flattening surfaces to 2D
- `Poisson Problem`: Solving PDEs on surfaces
- `Geodesic Distance`: Computing distances
- And more...

### Utils
Helper functions:
- `Solvers`: Eigenvalue problems and linear systems
- `MeshIO`: Reading/writing mesh files
- `Colormap`, `Distortion`: Visualization utilities

## Relationship to JavaScript Implementation

This Rust port maintains API compatibility with the original JavaScript implementation while leveraging Rust's performance and safety features. The documentation references the original JS code to help users familiar with the JavaScript version.

## Contributing

When adding new features:
1. Update the corresponding `.md` file in `docs/`
2. Follow the existing format: small code blocks with explanatory text
3. Include mathematical background where appropriate
4. Cross-reference related modules

## Additional Resources

- [Original JavaScript Framework](https://github.com/GeometryCollective/geometry-processing-js)
- [Discrete Differential Geometry Course](http://ddg.cs.cmu.edu/)
