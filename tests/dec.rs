use geometry_processing_rs::core::mesh::{Mesh};
use geometry_processing_rs::core::geometry::Geometry;
use geometry_processing_rs::core::dec::DEC;
use geometry_processing_rs::linear_algebra::{DenseMatrix};

mod common;
use common::{load_solution, parse_polygon_soup};

struct DecSolution {
    hodge0: Vec<f64>,
    hodge1: Vec<f64>,
    hodge2: Vec<f64>,
    phi: Vec<f64>,
    d_phi: Vec<f64>,
    omega: Vec<f64>,
    d_omega: Vec<f64>,
}

fn parse_dec_solution(content: &str) -> DecSolution {
    let mut hodge0 = Vec::new();
    let mut hodge1 = Vec::new();
    let mut hodge2 = Vec::new();
    let mut phi = Vec::new();
    let mut d_phi = Vec::new();
    let mut omega = Vec::new();
    let mut d_omega = Vec::new();

    for line in content.lines() {
        let line = line.trim();
        if line.is_empty() { continue; }
        let tokens: Vec<&str> = line.split_whitespace().collect();
        if tokens.len() < 2 { continue; }

        match tokens[0] {
            "hodge0" => hodge0.push(tokens[1].parse().unwrap()),
            "hodge1" => hodge1.push(tokens[1].parse().unwrap()),
            "hodge2" => hodge2.push(tokens[1].parse().unwrap()),
            "phi" => phi.push(tokens[1].parse().unwrap()),
            "dPhi" => d_phi.push(tokens[1].parse().unwrap()),
            "omega" => omega.push(tokens[1].parse().unwrap()),
            "dOmega" => d_omega.push(tokens[1].parse().unwrap()),
            _ => {}
        }
    }

    DecSolution { hodge0, hodge1, hodge2, phi, d_phi, omega, d_omega }
}

#[test]
fn test_dec() {
    let solution_data = load_solution("discrete-exterior-calculus");
    let soup = parse_polygon_soup(&solution_data);
    let mut mesh = Mesh::new();
    mesh.build(&soup);
    let geometry = Geometry::new(&mesh, soup.v, false);
    let sol = parse_dec_solution(&solution_data);

    let v_count = mesh.vertices.len();
    let e_count = mesh.edges.len();
    let f_count = mesh.faces.len();

    // HodgeStar0Form
    let h_star0 = DEC::build_hodge_star_0_form(&geometry);
    let h_star0_res = &h_star0 * &DenseMatrix::from_fn(v_count, 1, |_, _| 1.0);
    for i in 0..v_count {
        assert!((h_star0_res[(i, 0)] - sol.hodge0[i]).abs() < 1e-6);
    }

    // HodgeStar1Form
    let h_star1 = DEC::build_hodge_star_1_form(&geometry);
    let h_star1_res = &h_star1 * &DenseMatrix::from_fn(e_count, 1, |_, _| 1.0);
    for i in 0..e_count {
        assert!((h_star1_res[(i, 0)] - sol.hodge1[i]).abs() < 1e-6);
    }

    // HodgeStar2Form
    let h_star2 = DEC::build_hodge_star_2_form(&geometry);
    let h_star2_res = &h_star2 * &DenseMatrix::from_fn(f_count, 1, |_, _| 1.0);
    for i in 0..f_count {
        assert!((h_star2_res[(i, 0)] - sol.hodge2[i]).abs() < 1e-3);
    }

    // ExteriorDerivative0Form
    let d0 = DEC::build_exterior_derivative_0_form(&geometry);
    let mut phi_vec = DenseMatrix::zeros(v_count, 1);
    for i in 0..v_count { phi_vec[(i, 0)] = sol.phi[i]; }
    let d_phi_res = &d0 * &phi_vec;
    for i in 0..e_count {
        let diff = d_phi_res[(i, 0)] - sol.d_phi[i];
        let diff_neg = d_phi_res[(i, 0)] + sol.d_phi[i];
        assert!(diff.abs() < 1e-6 || diff_neg.abs() < 1e-6);
    }

    // ExteriorDerivative1Form
    let d1 = DEC::build_exterior_derivative_1_form(&geometry);
    let mut omega_vec = DenseMatrix::zeros(e_count, 1);
    for i in 0..e_count { omega_vec[(i, 0)] = sol.omega[i]; }
    let d_omega_res = &d1 * &omega_vec;
    for i in 0..f_count {
        let diff = d_omega_res[(i, 0)] - sol.d_omega[i];
        let diff_neg = d_omega_res[(i, 0)] + sol.d_omega[i];
        assert!(diff.abs() < 1e-6 || diff_neg.abs() < 1e-6);
    }

    // d1 * d0 = 0
    let d1d0 = &d1 * &d0;
    use geometry_processing_rs::linear_algebra::sparse_matrix::SparseMatrixMethods;
    assert!(d1d0.frobenius_norm() < 1e-6);
}
