use crate::core::mesh::{Mesh, MeshBackend};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Corner {
    pub index: usize,
}

impl Corner {
    pub fn new(index: usize) -> Self {
        Self { index }
    }

    pub fn vertex<B: MeshBackend>(&self, mesh: &Mesh<B>) -> usize {
        let h_idx = mesh.backend.corner_halfedge(self.index).expect("Corner should have a halfedge");
        let prev_idx = mesh.backend.halfedge_prev(h_idx).expect("Prev should exist");
        mesh.backend.halfedge_vertex(prev_idx).expect("Vertex should exist")
    }

    pub fn face<B: MeshBackend>(&self, mesh: &Mesh<B>) -> usize {
        let h_idx = mesh.backend.corner_halfedge(self.index).expect("Corner should have a halfedge");
        mesh.backend.halfedge_face(h_idx).expect("Face should exist")
    }

    pub fn next<B: MeshBackend>(&self, mesh: &Mesh<B>) -> Option<usize> {
        let h_idx = mesh.backend.corner_halfedge(self.index).expect("Corner should have a halfedge");
        let next_idx = mesh.backend.halfedge_next(h_idx).expect("Next should exist");
        mesh.backend.halfedge_corner(next_idx)
    }

    pub fn prev<B: MeshBackend>(&self, mesh: &Mesh<B>) -> Option<usize> {
        let h_idx = mesh.backend.corner_halfedge(self.index).expect("Corner should have a halfedge");
        let prev_idx = mesh.backend.halfedge_prev(h_idx).expect("Prev should exist");
        mesh.backend.halfedge_corner(prev_idx)
    }
}
