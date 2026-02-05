# Parameterization.rs Mapping

| Rust Function | JS Function | Description | Notes |
|---|---|---|---|
| `SpectralConformalParameterization::new` | (Constructor) | Initializes the spectral parameterizer. | |
| `SpectralConformalParameterization::build_conformal_energy` | `buildConformalEnergy` | Builds the complex conformal energy matrix. | |
| `SpectralConformalParameterization::flatten` | `flatten` | Computes the flattening via eigenvector of smallest eigenvalue. | |
| `BoundaryFirstFlattening::new` | (Constructor) | Initializes BFF. | |
| `BoundaryFirstFlattening::build_special_laplace` | | Builds the modified Laplacian for BFF. | |
| `BoundaryFirstFlattening::flatten` | `flatten` | Computes the flattening given boundary conditions. | |

# GeometricFlow.rs Mapping

| Rust Function | JS Function | Description | Notes |
|---|---|---|---|
| `MeanCurvatureFlow::integrate` | `integrate` | Performs one time step of mean curvature flow. | |
| `ModifiedMeanCurvatureFlow::integrate` | `integrate` | Performs MCF with fixed Laplacian structure. | |

# GeodesicDistance.rs Mapping
| Rust Function | JS Function | Description | Notes |
|---|---|---|---|
| `GeodesicDistance::compute` | `compute` | Computes geodesic distances using the Heat Method. | |

# Other Projects (Summary)
| File | Likely JS Counterpart | Description |
|---|---|---|
| `direction_field_design.rs` | `DirectionFieldDesign` | Algorithms for designing smooth direction fields. |
| `discrete_curvatures.rs` | `DiscreteCurvatures` | Visualization of discrete curvature measures. |
| `discrete_exterior_calculus_project.rs` | `DEC` project | Tests/Visualization for DEC operators. |
| `harmonic_bases.rs` | `HarmonicBases` | Construction of harmonic bases on surfaces. |
| `poisson_problem.rs` | `PoissonProblem` | Solving Poisson equation on meshes. |
| `simplicial_complex_operators.rs` | `SimplicialComplexOperators` | Topological operators (boundary, star, closure). |
| `tree_cotree.rs` | `TreeCotree` | Tree-Cotree decomposition for homology. |
| `vector_field_decomposition.rs` | `VectorFieldDecomposition` | Hodge-Helmholtz decomposition. |
