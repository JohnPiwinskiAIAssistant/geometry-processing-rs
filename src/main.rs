use bevy::prelude::*;
use bevy_mesh::{VertexAttributeValues, Indices, PrimitiveTopology};
use bevy::asset::RenderAssetUsages;
use geometry_processing_rs::core::mesh::{Mesh, PolygonSoup};
use geometry_processing_rs::core::geometry::Geometry;
use geometry_processing_rs::core::bevy_backend::BevyBackend;
use geometry_processing_rs::projects::geometric_flow::MeanCurvatureFlow;
use geometry_processing_rs::linear_algebra::Vector;
use std::fs;
use std::path::PathBuf;

#[derive(Resource)]
struct McfState {
    backend: BevyBackend,
    canonical_initial_positions: Vec<Vector>,
    canonical_positions: Vec<Vector>,
    flow: MeanCurvatureFlow,
    h: f64,
    mesh_handle: Handle<bevy::prelude::Mesh>,
    step_count: usize,
    is_running: bool,
    use_test_data: bool,
}

fn main() {
    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: "Mean Curvature Flow - Automated".to_string(),
                resolution: (1280, 720).into(),
                ..default()
            }),
            ..default()
        }))
        .add_systems(Startup, setup)
        .add_systems(Update, (run_mcf_system, sync_mesh_system))
        .run();
}

fn load_test_soup() -> PolygonSoup {
    let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    path.push("test-data");
    path.push("geometric-flow");
    path.push("solution.js");

    let content = fs::read_to_string(path).unwrap_or_else(|_| "".to_string());
    if content.is_empty() {
        return PolygonSoup { v: Vec::new(), f: Vec::new() };
    }
    
    let start = content.find('`').map(|i| i + 1).unwrap_or(0);
    let end = content.rfind('`').unwrap_or(content.len());
    let solution = &content[start..end];

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
                    positions.push(faer::mat![[x], [y], [z]]);
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

fn soup_to_bevy(soup: &PolygonSoup) -> bevy::prelude::Mesh {
    // In Bevy 0.14+ (0.18), we can use Default for RenderAssetUsages
    let mut mesh = bevy::prelude::Mesh::new(
        PrimitiveTopology::TriangleList,
        RenderAssetUsages::default(),
    );

    let positions: Vec<[f32; 3]> = soup.v.iter().map(|v| {
        [v[(0, 0)] as f32, v[(1, 0)] as f32, v[(2, 0)] as f32]
    }).collect();

    let indices = Indices::U32(soup.f.iter().map(|&i| i as u32).collect());

    mesh.insert_attribute(bevy::prelude::Mesh::ATTRIBUTE_POSITION, positions);
    mesh.insert_indices(indices);
    mesh.compute_flat_normals();
    mesh
}

fn setup(
    mut commands: Commands,
    mut meshes: ResMut<Assets<bevy::prelude::Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
) {
    // Start with Bevy primitive torus
    let torus = bevy::math::primitives::Torus {
        minor_radius: 0.3,
        major_radius: 1.0,
    };
    let bevy_mesh = bevy::prelude::Mesh::from(torus);
    let mesh_handle = meshes.add(bevy_mesh.clone());

    let backend = BevyBackend::new(bevy_mesh);
    let flow = MeanCurvatureFlow::new();

    let raw_pos_attr = backend.mesh.attribute(bevy::prelude::Mesh::ATTRIBUTE_POSITION).unwrap();
    let canonical_positions: Vec<Vector> = backend.canonical_to_raw.iter().map(|&raw_idx| {
        if let VertexAttributeValues::Float32x3(v) = raw_pos_attr {
            let p = v[raw_idx];
            faer::Mat::<f64>::from_fn(3, 1, |i, _| p[i] as f64)
        } else {
            panic!("Expected Float32x3 positions");
        }
    }).collect();
    
    println!("=== Mean Curvature Flow Demo ===");
    println!("Controls:");
    println!("  SPACE: Toggle simulation (Pause/Resume)");
    println!("  R    : Reset mesh");
    println!("  T    : Toggle Bevy Primitive / Test Solution Data");

    commands.insert_resource(McfState {
        backend,
        canonical_initial_positions: canonical_positions.clone(),
        canonical_positions,
        flow,
        h: 0.0001,
        mesh_handle: mesh_handle.clone(),
        step_count: 0,
        is_running: true,
        use_test_data: false,
    });

    commands.spawn((
        Mesh3d(mesh_handle),
        MeshMaterial3d(materials.add(StandardMaterial {
            base_color: Color::linear_rgb(0.3, 0.5, 0.9),
            metallic: 0.5,
            perceptual_roughness: 0.2,
            ..default()
        })),
        Transform::from_xyz(0.0, 0.0, 0.0),
    ));

    commands.spawn((
        PointLight {
            intensity: 1500.0,
            shadows_enabled: true,
            ..default()
        },
        Transform::from_xyz(4.0, 8.0, 4.0),
    ));

    commands.spawn((
        Camera3d::default(),
        Transform::from_xyz(0.0, 2.0, 5.0).looking_at(Vec3::ZERO, Vec3::Y),
    ));
}

fn run_mcf_system(
    mut mcf: ResMut<McfState>,
    mut meshes: ResMut<Assets<bevy::prelude::Mesh>>,
    keyboard: Res<ButtonInput<KeyCode>>,
) {
    if keyboard.just_pressed(KeyCode::Space) {
        mcf.is_running = !mcf.is_running;
        println!("MCF Simulation: {}", if mcf.is_running { "Running" } else { "Paused" });
    }

    if keyboard.just_pressed(KeyCode::KeyR) {
        mcf.canonical_positions = mcf.canonical_initial_positions.clone();
        mcf.step_count = 0;
        println!("MCF Simulation: Reset");
    }

    if keyboard.just_pressed(KeyCode::KeyT) {
        mcf.use_test_data = !mcf.use_test_data;
        let new_mesh = if mcf.use_test_data {
            println!("Switched to: Test Solution Torus");
            soup_to_bevy(&load_test_soup())
        } else {
            println!("Switched to: Bevy Primitive Torus");
            bevy::prelude::Mesh::from(bevy::math::primitives::Torus {
                minor_radius: 0.3,
                major_radius: 1.0,
            })
        };
        
        mcf.backend = BevyBackend::new(new_mesh.clone());
        let raw_pos_attr = mcf.backend.mesh.attribute(bevy::prelude::Mesh::ATTRIBUTE_POSITION).unwrap();
        mcf.canonical_positions = mcf.backend.canonical_to_raw.iter().map(|&raw_idx| {
            if let VertexAttributeValues::Float32x3(v) = raw_pos_attr {
                let p = v[raw_idx];
                faer::Mat::<f64>::from_fn(3, 1, |i, _| p[i] as f64)
            } else {
                panic!("Expected Float32x3 positions");
            }
        }).collect();
        mcf.canonical_initial_positions = mcf.canonical_positions.clone();
        mcf.step_count = 0;
        
        if let Some(mesh) = meshes.get_mut(&mcf.mesh_handle) {
            *mesh = new_mesh;
        }
    }

    if mcf.is_running {
        let mesh_wrapper = Mesh { backend: &mcf.backend };
        let mut geometry = Geometry {
            mesh: &mesh_wrapper,
            positions: mcf.canonical_positions.clone(),
        };

        let l = geometry.laplace_matrix();
        let m = geometry.mass_matrix();
        let lhs = &m - &(&l * mcf.h);
        
        if lhs.as_ref().sp_cholesky(faer::Side::Lower).is_ok() {
             mcf.flow.integrate(&mut geometry, mcf.h);
             mcf.canonical_positions = geometry.positions;
             mcf.step_count += 1;
             if mcf.step_count % 100 == 0 {
                 println!("MCF step: {}", mcf.step_count);
             }
        } else {
             println!("MCF Paused at step {}: Numerical instability detected.", mcf.step_count);
             mcf.is_running = false;
        }
    }
}

fn sync_mesh_system(
    mcf: Res<McfState>,
    mut meshes: ResMut<Assets<bevy::prelude::Mesh>>,
) {
    if let Some(mesh) = meshes.get_mut(&mcf.mesh_handle) {
        let mut new_positions = vec![[0.0f32; 3]; mcf.backend.mesh.count_vertices()];
        for (raw_idx, &canonical_idx) in mcf.backend.raw_to_canonical.iter().enumerate() {
            let p = &mcf.canonical_positions[canonical_idx];
            new_positions[raw_idx] = [p[(0, 0)] as f32, p[(1, 0)] as f32, p[(2, 0)] as f32];
        }

        mesh.insert_attribute(bevy::prelude::Mesh::ATTRIBUTE_POSITION, new_positions);
    }
}
