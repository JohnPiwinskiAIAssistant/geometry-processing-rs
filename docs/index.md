# Geometry Processing Rust Port Documentation

Welcome to the documentation for the Rust port of the `geometry-processing-js` framework. This port aims to provide a fast and flexible framework for 3D geometry processing, maintaining the same architectural principles as the original JavaScript implementation.

## Modules

### [Linear Algebra](./linear_algebra/index.md)
A wrapper around the `faer` and `num-complex` crates, providing a familiar API for complex numbers, vectors, and sparse/dense matrices.

### [Core](./core/index.md)
The heart of the framework, containing the halfedge mesh data structure and basic geometric operations.

### [Projects](./projects/index.md)
Implementations of various geometry processing algorithms, such as parameterization and Poisson problem solving.

### [Utils](./utils/index.md)
Common utilities for colormaps, mesh I/O, and solver wrappers.
