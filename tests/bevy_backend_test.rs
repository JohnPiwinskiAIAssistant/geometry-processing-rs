use bevy_mesh::{Mesh as BevyMesh, Indices, VertexAttributeValues};
use bevy::math::primitives::Torus;
use geometry_processing_rs::core::mesh::{Mesh, MeshBackend};
use geometry_processing_rs::core::geometry::Geometry;
use geometry_processing_rs::core::bevy_backend::BevyBackend;
use geometry_processing_rs::projects::geometric_flow::MeanCurvatureFlow;
use geometry_processing_rs::linear_algebra::Vector;

#[test]
fn test_bevy_mcf_headless() {
    // 1. Setup Bevy Mesh (Torus)
    let torus = Torus {
        minor_radius: 0.3,
        major_radius: 1.0,
    };
    let bevy_mesh = BevyMesh::from(torus);

    // 2. Extract positions
    let positions_attr = bevy_mesh.attribute(BevyMesh::ATTRIBUTE_POSITION).unwrap();
    let positions_vec = if let VertexAttributeValues::Float32x3(v) = positions_attr {
        v.iter().map(|p| {
            faer::Mat::<f64>::from_fn(3, 1, |i, _| p[i] as f64)
        }).collect::<Vec<Vector>>()
    } else {
        panic!("Position attribute not found or incorrect format");
    };

    // 3. Setup Backend and Flow
    println!("Building BevyBackend...");
    let backend = BevyBackend::new(bevy_mesh);
    let flow = MeanCurvatureFlow::new();

    let mesh_wrapper = Mesh { backend: &backend };
    let mut geometry = Geometry {
        mesh: &mesh_wrapper,
        positions: positions_vec,
    };

    // 4. Run Step (This is where the user thinks it fails)
    println!("Running MCF step...");
    flow.integrate(&mut geometry, 0.01);
    
    println!("Step successful. New position of vertex 0: {:?}", geometry.positions[0]);
}
