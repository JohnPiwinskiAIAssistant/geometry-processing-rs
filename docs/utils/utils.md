# Solvers.rs Mapping

| Rust Function | JS Function | JS Line Number | Description | Notes |
|---|---|---|---|---|
| `Solvers::residual` | (Various) | | Computes residual \|\|Ax - b\|\|. | Helper for validation |
| `Solvers::solve_inverse_power_method` | (Various) | | Solves generalized eigenvalue problem using inverse power method. | **Warning**: Uses dense matrix conversion (inefficient). |
| `Solvers::invert_2x2` | (Various) | | Inverts a 2x2 dense matrix. | Helper |

# Colormap.rs Mapping

| Rust Function | JS Function | Description | Notes |
|---|---|---|---|
| `Colormap::colormap` | `colormap` | Maps a value to a color using a seismic or other scheme. | |
| `Colormap::seismic_data` | | Returns data for seismic colormap. | |

# Distortion.rs Mapping

| Rust Function | JS Function | Description | Notes |
|---|---|---|---|
| `Distortion::compute_quasi_conformal_error_per_face` | | Computes quasi-conformal error (distortion) for a single face. | |
| `Distortion::compute_quasi_conformal_error` | | Computes total/average QC error and assigns colors to vertices. | |
| `Distortion::compute_area_scaling_per_face` | | Computes area scaling factor (log ratio). | |
| `hsv` | | Converts HSV to RGB Vector. | Helper |

# MeshIO.rs Mapping

| Rust Function | JS Function | Description | Notes |
|---|---|---|---|
| `MeshIO::read_obj` | `readOBJ` | Parses OBJ string into PolygonSoup. | Basic implementation, triangulated meshes only. |
| `MeshIO::write_obj` | | Serializes PolygonSoup to OBJ string. | |
