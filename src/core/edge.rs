use crate::core::mesh::Mesh;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Edge {
    pub halfedge: Option<usize>,
    pub index: usize,
}

impl Edge {
    pub fn new(index: usize) -> Self {
        Self {
            halfedge: None,
            index,
        }
    }

    pub fn on_boundary(&self, mesh: &Mesh) -> bool {
        if let Some(h_idx) = self.halfedge {
            let twin_idx = mesh.halfedges[h_idx].twin.expect("Twin should exist");
            mesh.halfedges[h_idx].on_boundary || mesh.halfedges[twin_idx].on_boundary
        } else {
            false
        }
    }
}
