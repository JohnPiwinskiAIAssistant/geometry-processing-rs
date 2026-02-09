use crate::linear_algebra::PolygonSoup;

pub struct MeshIO;

impl MeshIO {
    pub fn read_obj(input: &str) -> Option<PolygonSoup> {
        let mut positions = Vec::new();
        let mut indices = Vec::new();

        for line in input.lines() {
            let line = line.trim();
            if line.is_empty() { continue; }
            let tokens: Vec<&str> = line.split_whitespace().collect();
            if tokens.is_empty() { continue; }
            
            let identifier = tokens[0];

            match identifier {
                "v" => {
                    if tokens.len() >= 4 {
                        let x = tokens[1].parse::<f64>().ok()?;
                        let y = tokens[2].parse::<f64>().ok()?;
                        let z = tokens[3].parse::<f64>().ok()?;
                        positions.push(faer::mat![[x], [y], [z]]);
                    }
                }
                "f" => {
                    if tokens.len() > 4 {
                        eprintln!("Only triangle meshes are supported at this time!");
                        return None;
                    }
                    for i in 1..tokens.len() {
                        let index_str = tokens[i].split('/').next()?;
                        let index = index_str.parse::<usize>().ok()?;
                        indices.push(index - 1);
                    }
                }
                _ => {}
            }
        }

        Some(PolygonSoup { v: positions, f: indices })
    }

    pub fn write_obj(polygon_soup: &PolygonSoup) -> String {
        let mut output = String::new();

        // Write positions
        for p in &polygon_soup.v {
            output.push_str(&format!("v {} {} {}\n", p[(0, 0)], p[(1, 0)], p[(2, 0)]));
        }

        // Write indices
        for i in 0..(polygon_soup.f.len() / 3) {
            output.push_str("f ");
            for j in 0..3 {
                let index = polygon_soup.f[3 * i + j] + 1;
                output.push_str(&format!("{} ", index));
            }
            output.push('\n');
        }

        output
    }
}
