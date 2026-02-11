use crate::core::mesh::{Mesh, MeshBackend};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Face {
    pub index: usize,
}

impl Face {
    pub fn new(index: usize) -> Self {
        Self { index }
    }

    pub fn is_boundary_loop<B: MeshBackend>(&self, mesh: &Mesh<B>) -> bool {
        if let Some(h_idx) = mesh.backend.face_halfedge(self.index) {
            mesh.backend.halfedge_on_boundary(h_idx)
        } else {
            false
        }
    }

    pub fn adjacent_vertices<'a, B: MeshBackend>(&self, mesh: &'a Mesh<B>, ccw: bool) -> Box<dyn Iterator<Item = usize> + 'a> {
        mesh.face_adjacent_vertices(self.index, ccw)
    }

    pub fn adjacent_edges<'a, B: MeshBackend>(&self, mesh: &'a Mesh<B>, ccw: bool) -> Box<dyn Iterator<Item = usize> + 'a> {
        mesh.face_adjacent_edges(self.index, ccw)
    }

    pub fn adjacent_faces<'a, B: MeshBackend>(&self, mesh: &'a Mesh<B>, ccw: bool) -> Box<dyn Iterator<Item = usize> + 'a> {
        mesh.face_adjacent_faces(self.index, ccw)
    }

    pub fn adjacent_halfedges<'a, B: MeshBackend>(&self, mesh: &'a Mesh<B>, ccw: bool) -> Box<dyn Iterator<Item = usize> + 'a> {
        mesh.face_adjacent_halfedges(self.index, ccw)
    }

    pub fn adjacent_corners<'a, B: MeshBackend>(&self, mesh: &'a Mesh<B>, ccw: bool) -> Box<dyn Iterator<Item = usize> + 'a> {
        mesh.face_adjacent_corners(self.index, ccw)
    }
}
