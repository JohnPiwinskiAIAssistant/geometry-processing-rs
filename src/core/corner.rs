use crate::core::mesh::Mesh;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Corner {
    pub halfedge: Option<usize>,
    pub index: usize,
}

impl Corner {
    pub fn new(index: usize) -> Self {
        Self {
            halfedge: None,
            index,
        }
    }

    pub fn vertex(&self, mesh: &Mesh) -> usize {
        let h_idx = self.halfedge.expect("Corner should have a halfedge");
        let prev_idx = mesh.halfedges[h_idx].prev.expect("Prev should exist");
        mesh.halfedges[prev_idx].vertex.expect("Vertex should exist")
    }

    pub fn face(&self, mesh: &Mesh) -> usize {
        let h_idx = self.halfedge.expect("Corner should have a halfedge");
        mesh.halfedges[h_idx].face.expect("Face should exist")
    }

    pub fn next(&self, mesh: &Mesh) -> Option<usize> {
        let h_idx = self.halfedge.expect("Corner should have a halfedge");
        let next_idx = mesh.halfedges[h_idx].next.expect("Next should exist");
        mesh.halfedges[next_idx].corner
    }

    pub fn prev(&self, mesh: &Mesh) -> Option<usize> {
        let h_idx = self.halfedge.expect("Corner should have a halfedge");
        let prev_idx = mesh.halfedges[h_idx].prev.expect("Prev should exist");
        mesh.halfedges[prev_idx].corner
    }
}
