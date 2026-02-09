use geometry_processing_rs::core::mesh::{Mesh};
use geometry_processing_rs::core::geometry::Geometry;
use geometry_processing_rs::projects::geodesic_distance::{HeatMethod};
use geometry_processing_rs::linear_algebra::{DenseMatrix, Vector, Cholesky};
use geometry_processing_rs::linear_algebra::traits::{LinearSolver, SparseOps};

mod common;
use common::{load_solution, parse_polygon_soup};

struct GeodesicSolution {
    delta: DenseMatrix,
    x_sol: Vec<Vector>,
    div_sol: DenseMatrix,
    phi_sol: DenseMatrix,
}

fn parse_geodesic_solution(content: &str, v_count: usize, f_count: usize) -> GeodesicSolution {
    let mut delta = DenseMatrix::zeros(v_count, 1);
    let mut x_sol = vec![faer::Mat::zeros(3, 1); f_count];
    let mut div_sol = DenseMatrix::zeros(v_count, 1);
    let mut phi_sol = DenseMatrix::zeros(v_count, 1);

    let mut v_idx = 0;
    let mut f_idx = 0;
    let mut div_idx = 0;
    let mut phi_idx = 0;

    for line in content.lines() {
        let line = line.trim();
        if line.is_empty() { continue; }
        let tokens: Vec<&str> = line.split_whitespace().collect();
        if tokens.len() < 2 { continue; }

        match tokens[0] {
            "delta" => {
                delta[(v_idx, 0)] = tokens[1].parse().unwrap();
                v_idx += 1;
            }
            "X" => {
                x_sol[f_idx] = faer::mat![
                    [tokens[1].parse::<f64>().unwrap()],
                    [tokens[2].parse::<f64>().unwrap()],
                    [tokens[3].parse::<f64>().unwrap()],
                ];
                f_idx += 1;
            }
            "div" => {
                // solution file has a sign error
                div_sol[(div_idx, 0)] = -tokens[1].parse::<f64>().unwrap();
                div_idx += 1;
            }
            "phi" => {
                phi_sol[(phi_idx, 0)] = tokens[1].parse().unwrap();
                phi_idx += 1;
            }
            _ => {}
        }
    }

    GeodesicSolution { delta, x_sol, div_sol, phi_sol }
}

#[test]
fn test_geodesic_distance() {
    let solution_data = load_solution("geodesic-distance");
    let soup = parse_polygon_soup(&solution_data);
    let mut mesh = Mesh::new();
    mesh.build(&soup);
    let geometry = Geometry::new(&mesh, soup.v, false);
    
    let sol = parse_geodesic_solution(&solution_data, mesh.vertices.len(), mesh.faces.len());

    let heat_method = HeatMethod::new(&geometry);
    
    // computeVectorField
    let llt = Cholesky::new(&heat_method.f);
    let u = llt.solve(&sol.delta);
    println!("u min: {}, u max: {}", (0..u.nrows()).map(|i| u[(i, 0)]).fold(f64::INFINITY, f64::min), (0..u.nrows()).map(|i| u[(i, 0)]).fold(f64::NEG_INFINITY, f64::max));
    let x = heat_method.compute_vector_field(&u);
    for i in 0..mesh.faces.len() {
        let diff_mat = &x[i] - &sol.x_sol[i];
        let diff: f64 = (diff_mat.transpose() * &diff_mat).read(0, 0).sqrt();
        if diff >= 1e-5 {
            println!("Face {}: x={:?}, expected={:?}", i, x[i], sol.x_sol[i]);
            let normal = geometry.face_normal(&mesh.faces[i]).unwrap();
            println!("Normal: {:?}", normal);
            let area = geometry.area(&mesh.faces[i]);
            println!("Area: {}", area);
        }
        assert!(diff < 1e-5, "Vector field mismatch at face {}: got {:?}, expected {:?}", i, x[i], sol.x_sol[i]);
    }

    // computeDivergence
    let div = heat_method.compute_divergence(&x);
    for i in 0..mesh.vertices.len() {
        assert!((div[(i, 0)] - sol.div_sol[(i, 0)]).abs() < 1e-5, "Divergence mismatch at vertex {}: got {}, expected {}", i, div[(i, 0)], sol.div_sol[(i, 0)]);
    }

    // compute
    let phi = heat_method.compute(&sol.delta);
    for i in 0..mesh.vertices.len() {
        assert!((phi[(i, 0)] - sol.phi_sol[(i, 0)]).abs() < 1e-5, "Geodesic distance mismatch at vertex {}: got {}, expected {}", i, phi[(i, 0)], sol.phi_sol[(i, 0)]);
    }
}
