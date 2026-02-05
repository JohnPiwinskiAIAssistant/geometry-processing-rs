# Geometry.rs Mapping

| Rust Function | JS Function | JS Line Number | Description | Notes |
|---|---|---|---|---|
| `Geometry::new` | `Geometry.constructor` | 13 | Initializes Geometry for a mesh with vertex positions. | |
| `Geometry::normalize` | `normalize` | 387 | Centers mesh at origin and rescales to unit radius. | JS is global function |
| `Geometry::vector` | `Geometry.vector` | 33 | Computes the vector along a halfedge. | |
| `Geometry::length` | `Geometry.length` | 44 | Computes the length of an edge. | |
| `Geometry::midpoint` | `Geometry.midpoint` | 54 | Computes the midpoint of an edge. | |
| `Geometry::mean_edge_length` | `Geometry.meanEdgeLength` | 66 | Computes the mean edge length of the mesh. | |
| `Geometry::area` | `Geometry.area` | 80 | Computes the area of a face. | |
| `Geometry::total_area` | `Geometry.totalArea` | 93 | Computes the total surface area of the mesh. | |
| `Geometry::face_normal` | `Geometry.faceNormal` | 107 | Computes the normal of a face. | |
| `Geometry::centroid` | `Geometry.centroid` | 120 | Computes the centroid of a face. | |
| `Geometry::circumcenter` | `Geometry.circumcenter` | 136 | Computes the circumcenter of a face. | |
| `Geometry::orthonormal_bases` | `Geometry.orthonormalBases` | 158 | Computes an orthonormal basis for a face. | |
| `Geometry::angle` | `Geometry.angle` | 170 | Computes the angle (in radians) at a corner. | |
| `Geometry::cotan` | `Geometry.cotan` | 181 | Computes the cotangent of the angle opposite to a halfedge. | |
| `Geometry::dihedral_angle` | `Geometry.dihedralAngle` | 195 | Computes the signed dihedral angle between adjacent faces. | |
| `Geometry::barycentric_dual_area` | `Geometry.barycentricDualArea` | 210 | Computes the barycentric dual area of a vertex. | |
| `Geometry::circumcentric_dual_area` | `Geometry.circumcentricDualArea` | 224 | Computes the circumcentric dual area of a vertex. | |
| `Geometry::vertex_normal_equally_weighted` | `Geometry.vertexNormalEquallyWeighted` | 240 | Computes vertex normal using equal weights. | |
| `Geometry::vertex_normal_area_weighted` | `Geometry.vertexNormalAreaWeighted` | 256 | Computes vertex normal using face area weights. | |
| `Geometry::vertex_normal_angle_weighted` | `Geometry.vertexNormalAngleWeighted` | 272 | Computes vertex normal using tip angle weights. | |
| `Geometry::vertex_normal_gauss_curvature` | `Geometry.vertexNormalGaussCurvature` | 288 | Computes vertex normal using Gauss curvature method. | |
| `Geometry::vertex_normal_mean_curvature` | `Geometry.vertexNormalMeanCurvature` | 303 | Computes vertex normal using mean curvature method. | |
| `Geometry::vertex_normal_sphere_inscribed` | `Geometry.vertexNormalSphereInscribed` | 316 | Computes vertex normal using inscribed sphere method. | |
| `Geometry::angle_defect` | `Geometry.angleDefect` | 332 | Computes the angle defect at a vertex. | |
| `Geometry::scalar_gauss_curvature` | `Geometry.scalarGaussCurvature` | 345 | Computes scalar Gauss curvature at a vertex. | |
| `Geometry::scalar_mean_curvature` | `Geometry.scalarMeanCurvature` | 354 | Computes scalar mean curvature at a vertex. | |
| `Geometry::total_angle_defect` | `Geometry.totalAngleDefect` | 366 | Computes total angle defect (2π * Euler characteristic). | |
| `Geometry::principal_curvatures` | `Geometry.principalCurvatures` | 379 | Computes min and max principal curvatures at a vertex. | |
| `Geometry::laplace_matrix` | `Geometry.laplaceMatrix` | 398 | Builds the Cotan-Laplace matrix. | |
| `Geometry::mass_matrix` | `Geometry.massMatrix` | 419 | Builds the mass matrix (barycentric dual areas). | |
| `Geometry::complex_laplace_matrix` | `Geometry.complexLaplaceMatrix` | 433 | Builds the complex Cotan-Laplace matrix. | |
