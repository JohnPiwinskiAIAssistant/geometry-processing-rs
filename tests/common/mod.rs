use std::fs;
use std::path::PathBuf;
use geometry_processing_rs::core::mesh::PolygonSoup;
use geometry_processing_rs::linear_algebra::Vector;

pub fn load_solution(project: &str) -> String {
    let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    path.push("test-data");
    path.push(project);
    path.push("solution.js");

    let content = fs::read_to_string(path).expect("Failed to read solution file");
    
    // JS solution files start with: let solution = `
    // and end with: `;
    // Some might have extra content, but the data is inside backticks.
    let start = content.find('`').expect("Could not find start of solution data") + 1;
    let end = content.rfind('`').expect("Could not find end of solution data");
    
    content[start..end].to_string()
}

pub fn parse_polygon_soup(solution: &str) -> PolygonSoup {
    let mut positions = Vec::new();
    let mut indices = Vec::new();

    for line in solution.lines() {
        let tokens: Vec<&str> = line.split_whitespace().collect();
        if tokens.is_empty() { continue; }
        
        match tokens[0] {
            "v" => {
                if tokens.len() >= 4 {
                    let x = tokens[1].parse::<f64>().unwrap();
                    let y = tokens[2].parse::<f64>().unwrap();
                    let z = tokens[3].parse::<f64>().unwrap();
                    positions.push(Vector::new(x, y, z));
                }
            }
            "f" => {
                for i in 1..tokens.len() {
                    let index_str = tokens[i].split('/').next().unwrap();
                    let index = index_str.parse::<usize>().unwrap();
                    indices.push(index - 1);
                }
            }
            _ => {}
        }
    }

    PolygonSoup { v: positions, f: indices }
}

pub fn assert_near(a: f64, b: f64, eps: f64) {
    assert!((a - b).abs() < eps, "Value {} not near {} (diff: {})", a, b, (a - b).abs());
}

pub fn assert_vector_near(a: Vector, b: Vector, eps: f64) {
    let diff = a.minus(b).norm();
    assert!(diff < eps, "Vector {:?} not near {:?} (diff: {})", a, b, diff);
}
