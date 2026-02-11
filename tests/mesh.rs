use std::fs;
use std::path::PathBuf;
use geometry_processing_rs::core::mesh::{Mesh, MeshBackend, Vertex, Face};

mod common;
use common::{parse_polygon_soup};

fn load_input_mesh(_name: &str) -> String {
    let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    path.push("test-data");
    path.push("geometry");
    path.push("solution.js");

    let content = fs::read_to_string(path).expect("Failed to read input mesh file");
    
    let start = content.find('`').expect("Could not find start of mesh data") + 1;
    let end = content.rfind('`').expect("Could not find end of mesh data");
    
    content[start..end].to_string()
}

#[test]
fn test_mesh_connectivity() {
    let face_mesh_data = load_input_mesh("face");
    let soup = parse_polygon_soup(&face_mesh_data);
    let mut mesh = Mesh::new();
    let success = mesh.build(&soup);
    assert!(success, "Mesh build failed");

    // Vertex connectivity
    for v_idx in 0..mesh.num_vertices() {
        // vertex_adjacent_vertices
        let mut h_idx = mesh.backend.vertex_halfedge(v_idx).unwrap();
        for vv_idx in mesh.vertex_adjacent_vertices(v_idx, true) {
            let twin_idx = mesh.backend.halfedge_twin(h_idx).unwrap();
            let neighbor_idx = mesh.backend.halfedge_vertex(twin_idx).unwrap();
            assert_eq!(neighbor_idx, vv_idx);
            h_idx = mesh.backend.halfedge_next(twin_idx).unwrap();
        }

        h_idx = mesh.backend.vertex_halfedge(v_idx).unwrap();
        for vv_idx in mesh.vertex_adjacent_vertices(v_idx, false) {
            let twin_idx = mesh.backend.halfedge_twin(h_idx).unwrap();
            let neighbor_idx = mesh.backend.halfedge_vertex(twin_idx).unwrap();
            assert_eq!(neighbor_idx, vv_idx);
            let prev_idx = mesh.backend.halfedge_prev(h_idx).unwrap();
            h_idx = mesh.backend.halfedge_twin(prev_idx).unwrap();
        }

        // vertex_adjacent_edges
        h_idx = mesh.backend.vertex_halfedge(v_idx).unwrap();
        for e_idx in mesh.vertex_adjacent_edges(v_idx, true) {
            assert_eq!(mesh.backend.halfedge_edge(h_idx).unwrap(), e_idx);
            let twin_idx = mesh.backend.halfedge_twin(h_idx).unwrap();
            h_idx = mesh.backend.halfedge_next(twin_idx).unwrap();
        }

        // vertex_adjacent_faces
        h_idx = mesh.backend.vertex_halfedge(v_idx).unwrap();
        for f_idx in mesh.vertex_adjacent_faces(v_idx, true) {
            while mesh.backend.halfedge_on_boundary(h_idx) {
                let twin_idx = mesh.backend.halfedge_twin(h_idx).unwrap();
                h_idx = mesh.backend.halfedge_next(twin_idx).unwrap();
            }
            assert_eq!(mesh.backend.halfedge_face(h_idx).unwrap(), f_idx);
            let twin_idx = mesh.backend.halfedge_twin(h_idx).unwrap();
            h_idx = mesh.backend.halfedge_next(twin_idx).unwrap();
        }

        // vertex_adjacent_halfedges
        h_idx = mesh.backend.vertex_halfedge(v_idx).unwrap();
        for hh_idx in mesh.vertex_adjacent_halfedges(v_idx, true) {
            assert_eq!(h_idx, hh_idx);
            let twin_idx = mesh.backend.halfedge_twin(h_idx).unwrap();
            h_idx = mesh.backend.halfedge_next(twin_idx).unwrap();
        }

        // vertex_adjacent_corners
        h_idx = mesh.backend.vertex_halfedge(v_idx).unwrap();
        for c_idx in mesh.vertex_adjacent_corners(v_idx, true) {
            while mesh.backend.halfedge_on_boundary(h_idx) {
                let twin_idx = mesh.backend.halfedge_twin(h_idx).unwrap();
                h_idx = mesh.backend.halfedge_next(twin_idx).unwrap();
            }
            let next_h_idx = mesh.backend.halfedge_next(h_idx).unwrap();
            assert_eq!(mesh.backend.halfedge_corner(next_h_idx).unwrap(), c_idx);
            let twin_idx = mesh.backend.halfedge_twin(h_idx).unwrap();
            h_idx = mesh.backend.halfedge_next(twin_idx).unwrap();
        }
    }

    // Edge consistency
    for e_idx in 0..mesh.num_edges() {
        let h_idx = mesh.backend.edge_halfedge(e_idx).unwrap();
        assert_eq!(mesh.backend.halfedge_edge(h_idx).unwrap(), e_idx);
    }

    // Face connectivity
    for f_idx in 0..mesh.num_faces() {
        // face_adjacent_vertices
        let mut h_idx = mesh.backend.face_halfedge(f_idx).unwrap();
        for v_idx in mesh.face_adjacent_vertices(f_idx, true) {
            assert_eq!(mesh.backend.halfedge_vertex(h_idx).unwrap(), v_idx);
            h_idx = mesh.backend.halfedge_next(h_idx).unwrap();
        }

        // face_adjacent_edges
        h_idx = mesh.backend.face_halfedge(f_idx).unwrap();
        for e_idx in mesh.face_adjacent_edges(f_idx, true) {
            assert_eq!(mesh.backend.halfedge_edge(h_idx).unwrap(), e_idx);
            h_idx = mesh.backend.halfedge_next(h_idx).unwrap();
        }

        // face_adjacent_faces
        h_idx = mesh.backend.face_halfedge(f_idx).unwrap();
        for ff_idx in mesh.face_adjacent_faces(f_idx, true) {
            while mesh.backend.halfedge_on_boundary(mesh.backend.halfedge_twin(h_idx).unwrap()) {
                 h_idx = mesh.backend.halfedge_next(h_idx).unwrap();
            }
            let twin_h_idx = mesh.backend.halfedge_twin(h_idx).unwrap();
            assert_eq!(mesh.backend.halfedge_face(twin_h_idx).unwrap(), ff_idx);
            h_idx = mesh.backend.halfedge_next(h_idx).unwrap();
        }
    }

    // Corner consistency
    for c_idx in 0..mesh.num_corners() {
        let h_idx = mesh.backend.corner_halfedge(c_idx).unwrap();
        assert_eq!(mesh.backend.halfedge_corner(h_idx).unwrap(), c_idx);
    }
}
