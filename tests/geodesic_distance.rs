use geometry_processing_rs::core::mesh::{Mesh};
use geometry_processing_rs::core::geometry::Geometry;
use geometry_processing_rs::projects::geodesic_distance::{HeatMethod};
use geometry_processing_rs::linear_algebra::{DenseMatrix, Vector};

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
    let mut x_sol = vec![Vector::new(0.0, 0.0, 0.0); f_count];
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
                delta.set(tokens[1].parse().unwrap(), v_idx, 0);
                v_idx += 1;
            }
            "X" => {
                x_sol[f_idx] = Vector::new(
                    tokens[1].parse().unwrap(),
                    tokens[2].parse().unwrap(),
                    tokens[3].parse().unwrap(),
                );
                f_idx += 1;
            }
            "div" => {
                // solution file has a sign error
                div_sol.set(-tokens[1].parse::<f64>().unwrap(), div_idx, 0);
                div_idx += 1;
            }
            "phi" => {
                phi_sol.set(tokens[1].parse().unwrap(), phi_idx, 0);
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
    let llt = heat_method.f.chol();
    let u = llt.solve_positive_definite(&sol.delta);
    let x = heat_method.compute_vector_field(&u);
    for i in 0..mesh.faces.len() {
        let diff = x[i].minus(sol.x_sol[i]);
        assert!(diff.norm() < 1e-5, "Vector field mismatch at face {}: got {:?}, expected {:?}", i, x[i], sol.x_sol[i]);
    }

    // computeDivergence
    let div = heat_method.compute_divergence(&x);
    for i in 0..mesh.vertices.len() {
        assert!((div.get(i, 0) - sol.div_sol.get(i, 0)).abs() < 1e-5, "Divergence mismatch at vertex {}: got {}, expected {}", i, div.get(i, 0), sol.div_sol.get(i, 0));
    }

    // compute
    let phi = heat_method.compute(&sol.delta);
    for i in 0..mesh.vertices.len() {
        assert!((phi.get(i, 0) - sol.phi_sol.get(i, 0)).abs() < 1e-5, "Geodesic distance mismatch at vertex {}: got {}, expected {}", i, phi.get(i, 0), sol.phi_sol.get(i, 0));
    }
}
