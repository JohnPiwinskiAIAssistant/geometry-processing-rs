use geometry_processing_rs::core::mesh::{Mesh};
use geometry_processing_rs::core::geometry::Geometry;
use geometry_processing_rs::projects::geometric_flow::{MeanCurvatureFlow, ModifiedMeanCurvatureFlow};
use geometry_processing_rs::linear_algebra::{Vector};

mod common;
use common::{load_solution, parse_polygon_soup};

struct FlowSolution {
    steps: usize,
    h: f64,
    mcf_positions: Vec<Vector>,
    mmcf_positions: Vec<Vector>,
}

fn parse_flow_solution(content: &str) -> FlowSolution {
    let mut steps = 0;
    let mut h = 0.0;
    let mut mcf_positions = Vec::new();
    let mut mmcf_positions = Vec::new();

    for line in content.lines() {
        let line = line.trim();
        if line.is_empty() { continue; }
        let tokens: Vec<&str> = line.split_whitespace().collect();
        if tokens.len() < 2 { continue; }

        match tokens[0] {
            "steps" => steps = tokens[1].parse().unwrap(),
            "h" => h = tokens[1].parse().unwrap(),
            "mcf" => {
                mcf_positions.push(Vector::new(
                    tokens[1].parse().unwrap(),
                    tokens[2].parse().unwrap(),
                    tokens[3].parse().unwrap(),
                ));
            }
            "mmcf" => {
                mmcf_positions.push(Vector::new(
                    tokens[1].parse().unwrap(),
                    tokens[2].parse().unwrap(),
                    tokens[3].parse().unwrap(),
                ));
            }
            _ => {}
        }
    }

    FlowSolution { steps, h, mcf_positions, mmcf_positions }
}

#[test]
fn test_geometric_flow() {
    let solution_data = load_solution("geometric-flow");
    let soup = parse_polygon_soup(&solution_data);
    let mut mesh = Mesh::new();
    mesh.build(&soup);
    
    let sol = parse_flow_solution(&solution_data);

    // Test MCF
    {
        let mut geometry = Geometry::new(&mesh, soup.v.clone(), false);
        {
            let mut flow = MeanCurvatureFlow::new(&mut geometry);
            for _ in 0..sol.steps {
                flow.integrate(sol.h);
            }
        } // flow dropped here

        for i in 0..sol.mcf_positions.len() {
            let diff = geometry.positions[i].minus(sol.mcf_positions[i]);
            assert!(diff.norm() < 1e-4, "MCF mismatch at vertex {}: got {:?}, expected {:?}", i, geometry.positions[i], sol.mcf_positions[i]);
        }
    }

    // Test MMCF
    {
        let mut geometry = Geometry::new(&mesh, soup.v, false);
        {
            let mut flow = ModifiedMeanCurvatureFlow::new(&mut geometry);
            for _ in 0..sol.steps {
                flow.integrate(sol.h);
            }
        } // flow dropped here

        for i in 0..sol.mmcf_positions.len() {
            let diff = geometry.positions[i].minus(sol.mmcf_positions[i]);
            assert!(diff.norm() < 1e-4, "MMCF mismatch at vertex {}: got {:?}, expected {:?}", i, geometry.positions[i], sol.mmcf_positions[i]);
        }
    }
}
