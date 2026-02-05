use crate::core::mesh::Mesh;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Face {
    pub halfedge: Option<usize>,
    pub index: usize,
}

impl Face {
    pub fn new(index: usize) -> Self {
        Self {
            halfedge: None,
            index,
        }
    }

    pub fn is_boundary_loop(&self, mesh: &Mesh) -> bool {
        if let Some(h_idx) = self.halfedge {
            mesh.halfedges[h_idx].on_boundary
        } else {
            false
        }
    }

    pub fn adjacent_vertices<'a>(&self, mesh: &'a Mesh, ccw: bool) -> impl Iterator<Item = usize> + 'a {
        mesh.face_adjacent_vertices(self.index, ccw)
    }

    pub fn adjacent_edges<'a>(&self, mesh: &'a Mesh, ccw: bool) -> impl Iterator<Item = usize> + 'a {
        mesh.face_adjacent_edges(self.index, ccw)
    }

    pub fn adjacent_faces<'a>(&self, mesh: &'a Mesh, ccw: bool) -> impl Iterator<Item = usize> + 'a {
        mesh.face_adjacent_faces(self.index, ccw)
    }

    pub fn adjacent_halfedges<'a>(&self, mesh: &'a Mesh, ccw: bool) -> impl Iterator<Item = usize> + 'a {
        mesh.face_adjacent_halfedges(self.index, ccw)
    }

    pub fn adjacent_corners<'a>(&self, mesh: &'a Mesh, ccw: bool) -> impl Iterator<Item = usize> + 'a {
        mesh.face_adjacent_corners(self.index, ccw)
    }
}
