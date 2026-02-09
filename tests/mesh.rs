use std::fs;
use std::path::PathBuf;
use geometry_processing_rs::core::mesh::Mesh;

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
    for v in &mesh.vertices {
        let v_idx = v.index;
        
        // vertex_adjacent_vertices
        let mut h_idx = mesh.vertices[v_idx].halfedge.unwrap();
        for vv_idx in mesh.vertex_adjacent_vertices(v_idx, true) {
            let twin_idx = mesh.halfedges[h_idx].twin.unwrap();
            let neighbor_idx = mesh.halfedges[twin_idx].vertex.unwrap();
            assert_eq!(neighbor_idx, vv_idx);
            h_idx = mesh.halfedges[twin_idx].next.unwrap();
        }

        h_idx = mesh.vertices[v_idx].halfedge.unwrap();
        for vv_idx in mesh.vertex_adjacent_vertices(v_idx, false) {
            let twin_idx = mesh.halfedges[h_idx].twin.unwrap();
            let neighbor_idx = mesh.halfedges[twin_idx].vertex.unwrap();
            assert_eq!(neighbor_idx, vv_idx);
            h_idx = mesh.halfedges[mesh.halfedges[h_idx].prev.unwrap()].twin.unwrap();
        }

        // vertex_adjacent_edges
        h_idx = mesh.vertices[v_idx].halfedge.unwrap();
        for e_idx in mesh.vertex_adjacent_edges(v_idx, true) {
            assert_eq!(mesh.halfedges[h_idx].edge.unwrap(), e_idx);
            let twin_idx = mesh.halfedges[h_idx].twin.unwrap();
            h_idx = mesh.halfedges[twin_idx].next.unwrap();
        }

        // vertex_adjacent_faces
        h_idx = mesh.vertices[v_idx].halfedge.unwrap();
        for f_idx in mesh.vertex_adjacent_faces(v_idx, true) {
            while mesh.halfedges[h_idx].on_boundary {
                let twin_idx = mesh.halfedges[h_idx].twin.unwrap();
                h_idx = mesh.halfedges[twin_idx].next.unwrap();
            }
            assert_eq!(mesh.halfedges[h_idx].face.unwrap(), f_idx);
            let twin_idx = mesh.halfedges[h_idx].twin.unwrap();
            h_idx = mesh.halfedges[twin_idx].next.unwrap();
        }

        // vertex_adjacent_halfedges
        h_idx = mesh.vertices[v_idx].halfedge.unwrap();
        for hh_idx in mesh.vertex_adjacent_halfedges(v_idx, true) {
            assert_eq!(h_idx, hh_idx);
            let twin_idx = mesh.halfedges[h_idx].twin.unwrap();
            h_idx = mesh.halfedges[twin_idx].next.unwrap();
        }

        // vertex_adjacent_corners
        h_idx = mesh.vertices[v_idx].halfedge.unwrap();
        for c_idx in mesh.vertex_adjacent_corners(v_idx, true) {
            while mesh.halfedges[h_idx].on_boundary {
                let twin_idx = mesh.halfedges[h_idx].twin.unwrap();
                h_idx = mesh.halfedges[twin_idx].next.unwrap();
            }
            let next_h_idx = mesh.halfedges[h_idx].next.unwrap();
            assert_eq!(mesh.halfedges[next_h_idx].corner.unwrap(), c_idx);
            let twin_idx = mesh.halfedges[h_idx].twin.unwrap();
            h_idx = mesh.halfedges[twin_idx].next.unwrap();
        }
    }

    // Edge consistency
    for e in &mesh.edges {
        let h_idx = e.halfedge.unwrap();
        assert_eq!(mesh.halfedges[h_idx].edge.unwrap(), e.index);
    }

    // Face connectivity
    for f in &mesh.faces {
        let f_idx = f.index;
        
        // face_adjacent_vertices
        let mut h_idx = mesh.faces[f_idx].halfedge.unwrap();
        for v_idx in mesh.face_adjacent_vertices(f_idx, true) {
            assert_eq!(mesh.halfedges[h_idx].vertex.unwrap(), v_idx);
            h_idx = mesh.halfedges[h_idx].next.unwrap();
        }

        // face_adjacent_edges
        h_idx = mesh.faces[f_idx].halfedge.unwrap();
        for e_idx in mesh.face_adjacent_edges(f_idx, true) {
            assert_eq!(mesh.halfedges[h_idx].edge.unwrap(), e_idx);
            h_idx = mesh.halfedges[h_idx].next.unwrap();
        }

        // face_adjacent_faces
        h_idx = mesh.faces[f_idx].halfedge.unwrap();
        for ff_idx in mesh.face_adjacent_faces(f_idx, true) {
            while mesh.halfedges[mesh.halfedges[h_idx].twin.unwrap()].on_boundary {
                 h_idx = mesh.halfedges[h_idx].next.unwrap();
            }
            let twin_h_idx = mesh.halfedges[h_idx].twin.unwrap();
            assert_eq!(mesh.halfedges[twin_h_idx].face.unwrap(), ff_idx);
            h_idx = mesh.halfedges[h_idx].next.unwrap();
        }
    }

    // Corner consistency
    for c in &mesh.corners {
        let h_idx = c.halfedge.unwrap();
        assert_eq!(mesh.halfedges[h_idx].corner.unwrap(), c.index);
    }
}
