use bevy_mesh::Mesh as BevyMesh;
use geometry_processing_rs::core::bevy_backend::BevyBackend;
use geometry_processing_rs::core::mesh::MeshBackend;

fn main() {
    println!("Checking Bevy torus topology...");
    let torus = bevy::math::primitives::Torus {
        minor_radius: 0.3,
        major_radius: 1.0,
    };
    let bevy_mesh = BevyMesh::from(torus);
    let backend = BevyBackend::new(bevy_mesh);
    
    let v = backend.num_vertices();
    let e = backend.num_edges();
    let f = backend.num_faces();
    let h = backend.num_halfedges();
    
    // Euler characteristic: V - E + F
    let chi = (v as i32) - (e as i32) + (f as i32);
    
    println!("Vertices (V): {}", v);
    println!("Edges (E):    {}", e);
    println!("Faces (F):    {}", f);
    println!("Halfedges:    {}", h);
    println!("Euler characteristic (χ = V - E + F): {}", chi);
    
    // Boundaries
    let mut boundary_count = 0;
    for h_idx in 0..h {
        if backend.halfedge_on_boundary(h_idx) {
            boundary_count += 1;
        }
    }
    println!("Boundary half-edges: {}", boundary_count);

    if chi == 0 && boundary_count == 0 {
        println!("\nVERIFIED: The mesh is a closed manifold with genus 1 (a true torus).");
    } else if boundary_count > 0 {
        println!("\nWARNING: The mesh has boundaries! It is an open sheet.");
    } else {
        println!("\nWARNING: Unexpected Euler characteristic: {}", chi);
    }
}
