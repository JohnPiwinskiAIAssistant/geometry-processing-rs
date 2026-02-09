use geometry_processing_rs::core::mesh::{Mesh};
use geometry_processing_rs::core::geometry::Geometry;
use geometry_processing_rs::projects::parameterization::{SpectralConformalParameterization};
use geometry_processing_rs::linear_algebra::{Complex, ComplexTriplet, ComplexSparseMatrix, Vector};
use num_complex::Complex64;

mod common;
use common::{load_solution, parse_polygon_soup};

struct ParameterizationSolution {
    ec_sol: ComplexSparseMatrix,
    uv_sol: Vec<Vector>,
}

fn parse_parameterization_solution(content: &str, v_count: usize) -> ParameterizationSolution {
    let mut t = ComplexTriplet::new(v_count, v_count);
    let mut uv_sol = vec![faer::Mat::zeros(3, 1); v_count];

    let mut uv_idx = 0;

    for line in content.lines() {
        let line = line.trim();
        if line.is_empty() { continue; }
        let tokens: Vec<&str> = line.split_whitespace().collect();
        if tokens.len() < 2 { continue; }

        match tokens[0] {
            "T" => {
                t.add_entry(
                    Complex::new(tokens[1].parse().unwrap(), tokens[2].parse().unwrap()),
                    tokens[3].parse().unwrap(),
                    tokens[4].parse().unwrap(),
                );
            }
            "uv" => {
                uv_sol[uv_idx] = faer::mat![
                    [tokens[1].parse::<f64>().unwrap()],
                    [tokens[2].parse::<f64>().unwrap()],
                    [0.0],
                ];
                uv_idx += 1;
            }
            _ => {}
        }
    }

    use geometry_processing_rs::linear_algebra::sparse_matrix::SparseMatrixMethods;
    ParameterizationSolution {
        ec_sol: ComplexSparseMatrix::from_triplets(v_count, v_count, &t.data),
        uv_sol,
    }
}

#[test]
fn test_parameterization() {
    let solution_data = load_solution("parameterization");
    let soup = parse_polygon_soup(&solution_data);
    println!("Soup has {} vertices, {} faces", soup.v.len(), soup.f.len());
    let mut mesh = Mesh::new();
    mesh.build(&soup);
    let geometry = Geometry::new(&mesh, soup.v, false);
    
    let sol = parse_parameterization_solution(&solution_data, mesh.vertices.len());

    let scp = SpectralConformalParameterization::new(&geometry);
    
    println!("Number of boundaries: {}", geometry.mesh.boundaries.len());
    for (i, b) in geometry.mesh.boundaries.iter().enumerate() {
        let count = geometry.mesh.face_adjacent_halfedges(b.index, true).count();
        println!("Boundary {} has {} edges", i, count);
    }
    let ec = scp.build_conformal_energy();
    
    use geometry_processing_rs::linear_algebra::sparse_matrix::SparseMatrixMethods;
    // Debug: check norms and NNZ
    println!("EC norm: {}, NNZ: {}", ec.frobenius_norm(), ec.compute_nnz());
    println!("Sol EC norm: {}, NNZ: {}", sol.ec_sol.frobenius_norm(), sol.ec_sol.compute_nnz());
    
    let mut ed = geometry.complex_laplace_matrix();
    ed = ed.scale(Complex64::new(0.5, 0.0));
    println!("My ED norm: {}", ed.frobenius_norm());
    
    // Extract real part from sol.ec_sol
    let mut sol_ed_triplet = ComplexTriplet::new(mesh.vertices.len(), mesh.vertices.len());
    let sol_values = sol.ec_sol.values();
    for j in 0..sol.ec_sol.ncols() {
        let start = sol.ec_sol.col_ptrs()[j];
        let end = sol.ec_sol.col_ptrs()[j+1];
        for k in start..end {
            sol_ed_triplet.add_entry(Complex64::new(sol_values.re[k], 0.0), sol.ec_sol.row_indices()[k], j);
        }
    }
    let sol_ed = ComplexSparseMatrix::from_triplets(mesh.vertices.len(), mesh.vertices.len(), &sol_ed_triplet.data);
    println!("Sol ED norm: {}", sol_ed.frobenius_norm());
    
    let a = &ec - &ed;
    println!("My A norm: {}", a.frobenius_norm());
    
    // Extract imag part from sol.ec_sol
    let mut sol_a_triplet = ComplexTriplet::new(mesh.vertices.len(), mesh.vertices.len());
    for j in 0..sol.ec_sol.ncols() {
        let start = sol.ec_sol.col_ptrs()[j];
        let end = sol.ec_sol.col_ptrs()[j+1];
        for k in start..end {
            sol_a_triplet.add_entry(Complex64::new(0.0, sol_values.im[k]), sol.ec_sol.row_indices()[k], j);
        }
    }
    let sol_a = ComplexSparseMatrix::from_triplets(mesh.vertices.len(), mesh.vertices.len(), &sol_a_triplet.data);
    println!("Sol A norm: {}", sol_a.frobenius_norm());
    
    let total_diff = (&sol.ec_sol - &ec).frobenius_norm();
    assert!(total_diff < 1e-6, "Conformal energy matrix mismatch: norm = {}", total_diff);

    // flatten
    let uv = scp.flatten();
    for i in 0..mesh.vertices.len() {
        let diff_mat = &uv[i] - &sol.uv_sol[i];
        let diff: f64 = (diff_mat.transpose() * &diff_mat).read(0, 0).sqrt();
        assert!(diff < 1.0, "Flattening mismatch at vertex {}: got {:?}, expected {:?}", i, uv[i], sol.uv_sol[i]);
    }
}
