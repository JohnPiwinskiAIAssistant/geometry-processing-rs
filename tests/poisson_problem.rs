use geometry_processing_rs::core::mesh::{Mesh};
use geometry_processing_rs::core::geometry::Geometry;
use geometry_processing_rs::projects::poisson_problem::{ScalarPoissonProblem};
use geometry_processing_rs::linear_algebra::{DenseMatrix};

mod common;
use common::{load_solution, parse_polygon_soup};

struct PoissonSolution {
    rho: DenseMatrix,
    phi_sol: DenseMatrix,
}

fn parse_poisson_solution(content: &str, v_count: usize) -> PoissonSolution {
    let mut rho = DenseMatrix::zeros(v_count, 1);
    let mut phi_sol = DenseMatrix::zeros(v_count, 1);

    let mut rho_idx = 0;
    let mut phi_idx = 0;

    for line in content.lines() {
        let line = line.trim();
        if line.is_empty() { continue; }
        let tokens: Vec<&str> = line.split_whitespace().collect();
        if tokens.len() < 2 { continue; }

        match tokens[0] {
            "rho" => {
                rho.set(tokens[1].parse().unwrap(), rho_idx, 0);
                rho_idx += 1;
            }
            "phi" => {
                phi_sol.set(tokens[1].parse().unwrap(), phi_idx, 0);
                phi_idx += 1;
            }
            _ => {}
        }
    }

    PoissonSolution { rho, phi_sol }
}

#[test]
fn test_poisson_problem() {
    let solution_data = load_solution("poisson-problem");
    let soup = parse_polygon_soup(&solution_data);
    let mut mesh = Mesh::new();
    mesh.build(&soup);
    let geometry = Geometry::new(&mesh, soup.v, false);
    
    let sol = parse_poisson_solution(&solution_data, mesh.vertices.len());

    let scalar_poisson_problem = ScalarPoissonProblem::new(&geometry);
    let phi = scalar_poisson_problem.solve(&sol.rho);

    for i in 0..mesh.vertices.len() {
        assert!((phi.get(i, 0) - sol.phi_sol.get(i, 0)).abs() < 1e-3, "Poisson solution mismatch at vertex {}: got {}, expected {}", i, phi.get(i, 0), sol.phi_sol.get(i, 0));
    }
}
