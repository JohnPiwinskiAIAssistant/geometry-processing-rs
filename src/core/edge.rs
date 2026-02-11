use crate::core::mesh::{Mesh, MeshBackend};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Edge {
    pub index: usize,
}

impl Edge {
    pub fn new(index: usize) -> Self {
        Self { index }
    }

    pub fn on_boundary<B: MeshBackend>(&self, mesh: &Mesh<B>) -> bool {
        if let Some(h_idx) = mesh.backend.edge_halfedge(self.index) {
            let twin_idx = mesh.backend.halfedge_twin(h_idx).expect("Twin should exist");
            mesh.backend.halfedge_on_boundary(h_idx) || mesh.backend.halfedge_on_boundary(twin_idx)
        } else {
            false
        }
    }
}
