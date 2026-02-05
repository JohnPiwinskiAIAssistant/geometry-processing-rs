use geometry_processing_rs::core::mesh::Mesh;
use geometry_processing_rs::core::geometry::Geometry;
use geometry_processing_rs::linear_algebra::{Vector, DenseMatrix};

mod common;
use common::{load_solution, parse_polygon_soup, assert_near, assert_vector_near};

#[test]
fn test_geometry() {
    let solution = load_solution("geometry");
    let soup = parse_polygon_soup(&solution);
    let mut mesh = Mesh::new();
    mesh.build(&soup);
    let geometry = Geometry::new(&mesh, soup.v, false);

    let mut barycentric_dual_areas_sol = Vec::new();
    let mut circumcentric_dual_areas_sol = Vec::new();
    let mut angle_weighted_normals_sol = Vec::new();
    let mut sphere_inscribed_normals_sol = Vec::new();
    let mut area_weighted_normals_sol = Vec::new();
    let mut gauss_curvature_normals_sol = Vec::new();
    let mut mean_curvature_normals_sol = Vec::new();
    let mut scalar_gauss_curvatures_sol = Vec::new();
    let mut scalar_mean_curvatures_sol = Vec::new();
    let mut k1_sol = Vec::new();
    let mut k2_sol = Vec::new();
    let mut total_angle_defect_sol = 0.0;
    
    let mut laplace_triplets = Vec::new();
    let mut mass_sol = Vec::new();

    for line in solution.lines() {
        let tokens: Vec<&str> = line.split_whitespace().collect();
        if tokens.is_empty() { continue; }
        
        match tokens[0] {
            "barycentricDualArea" => barycentric_dual_areas_sol.push(tokens[1].parse::<f64>().unwrap()),
            "circumcentricDualArea" => circumcentric_dual_areas_sol.push(tokens[1].parse::<f64>().unwrap()),
            "angleWeightedNormal" => angle_weighted_normals_sol.push(Vector::new(
                tokens[1].parse().unwrap(), tokens[2].parse().unwrap(), tokens[3].parse().unwrap()
            )),
            "sphereInscribedNormal" => sphere_inscribed_normals_sol.push(Vector::new(
                tokens[1].parse().unwrap(), tokens[2].parse().unwrap(), tokens[3].parse().unwrap()
            )),
            "areaWeightedNormal" => area_weighted_normals_sol.push(Vector::new(
                tokens[1].parse().unwrap(), tokens[2].parse().unwrap(), tokens[3].parse().unwrap()
            )),
            "gaussCurvatureNormal" => gauss_curvature_normals_sol.push(Vector::new(
                tokens[1].parse().unwrap(), tokens[2].parse().unwrap(), tokens[3].parse().unwrap()
            )),
            "meanCurvatureNormal" => mean_curvature_normals_sol.push(Vector::new(
                tokens[1].parse().unwrap(), tokens[2].parse().unwrap(), tokens[3].parse().unwrap()
            )),
            "scalarGaussCurvature" => scalar_gauss_curvatures_sol.push(tokens[1].parse::<f64>().unwrap()),
            "scalarMeanCurvature" => scalar_mean_curvatures_sol.push(tokens[1].parse::<f64>().unwrap()),
            "k1" => k1_sol.push(tokens[1].parse::<f64>().unwrap()),
            "k2" => k2_sol.push(tokens[1].parse::<f64>().unwrap()),
            "totalAngleDefect" => total_angle_defect_sol = tokens[1].parse::<f64>().unwrap(),
            "T" => laplace_triplets.push((
                tokens[1].parse::<f64>().unwrap(),
                tokens[2].parse::<usize>().unwrap(),
                tokens[3].parse::<usize>().unwrap()
            )),
            "mass" => mass_sol.push(tokens[1].parse::<f64>().unwrap()),
            _ => {}
        }
    }

    // Test Dual Areas
    for v in &mesh.vertices {
        assert_near(barycentric_dual_areas_sol[v.index], geometry.barycentric_dual_area(v), 1e-5);
        assert_near(circumcentric_dual_areas_sol[v.index], geometry.circumcentric_dual_area(v), 1e-5);
    }

    // Test Normals
    for v in &mesh.vertices {
        assert_vector_near(angle_weighted_normals_sol[v.index], geometry.vertex_normal_angle_weighted(v), 1e-5);
        assert_vector_near(sphere_inscribed_normals_sol[v.index], geometry.vertex_normal_sphere_inscribed(v), 1e-5);
        assert_vector_near(area_weighted_normals_sol[v.index], geometry.vertex_normal_area_weighted(v), 1e-5);
        
        let gn = geometry.vertex_normal_gauss_curvature(v);
        let sol_gn = gauss_curvature_normals_sol[v.index];
        if gn.minus(sol_gn).norm() > 1e-5 && gn.negated().minus(sol_gn).norm() > 1e-5 {
             panic!("Gauss normal mismatch at vertex {}", v.index);
        }

        let mn = geometry.vertex_normal_mean_curvature(v);
        let sol_mn = mean_curvature_normals_sol[v.index];
        if mn.minus(sol_mn).norm() > 1e-5 && mn.negated().minus(sol_mn).norm() > 1e-5 {
             panic!("Mean normal mismatch at vertex {}", v.index);
        }
    }

    // Test Curvatures
    for v in &mesh.vertices {
        assert_near(scalar_gauss_curvatures_sol[v.index], geometry.scalar_gauss_curvature(v), 1e-5);
        assert_near(scalar_mean_curvatures_sol[v.index], geometry.scalar_mean_curvature(v), 1e-5);
        
        let [mut k1, mut k2] = geometry.principal_curvatures(v);
        if k1.abs() > k2.abs() { std::mem::swap(&mut k1, &mut k2); }
        assert_near(k1_sol[v.index], k1, 1e-5);
        assert_near(k2_sol[v.index], k2, 1e-5);
    }

    assert_near(total_angle_defect_sol, geometry.total_angle_defect(), 1e-5);

    // Test Laplace Matrix
    let laplace = geometry.laplace_matrix();
    
    for (val, r, c) in laplace_triplets {
        assert_near(val, laplace.get(r, c), 1e-5);
    }

    // Test Mass Matrix
    let mass_matrix = geometry.mass_matrix();
    let ones = DenseMatrix::ones(mesh.vertices.len(), 1);
    let diag = mass_matrix.times_dense(&ones);
    for i in 0..mesh.vertices.len() {
        assert_near(mass_sol[i], diag.get(i, 0), 1e-5);
    }
}
