use crate::core::mesh::Mesh;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Vertex {
    pub halfedge: Option<usize>,
    pub index: usize,
}

impl Vertex {
    pub fn new(index: usize) -> Self {
        Self {
            halfedge: None,
            index,
        }
    }

    pub fn is_isolated(&self) -> bool {
        self.halfedge.is_none()
    }

    pub fn on_boundary(&self, mesh: &Mesh) -> bool {
        mesh.on_boundary(self.index)
    }

    pub fn degree(&self, mesh: &Mesh) -> usize {
        mesh.vertex_degree(self.index)
    }

    pub fn adjacent_vertices<'a>(&self, mesh: &'a Mesh, ccw: bool) -> impl Iterator<Item = usize> + 'a {
        mesh.vertex_adjacent_vertices(self.index, ccw)
    }

    pub fn adjacent_edges<'a>(&self, mesh: &'a Mesh, ccw: bool) -> impl Iterator<Item = usize> + 'a {
        mesh.vertex_adjacent_edges(self.index, ccw)
    }

    pub fn adjacent_faces<'a>(&self, mesh: &'a Mesh, ccw: bool) -> impl Iterator<Item = usize> + 'a {
        mesh.vertex_adjacent_faces(self.index, ccw)
    }

    pub fn adjacent_halfedges<'a>(&self, mesh: &'a Mesh, ccw: bool) -> impl Iterator<Item = usize> + 'a {
        mesh.vertex_adjacent_halfedges(self.index, ccw)
    }

    pub fn adjacent_corners<'a>(&self, mesh: &'a Mesh, ccw: bool) -> impl Iterator<Item = usize> + 'a {
        mesh.vertex_adjacent_corners(self.index, ccw)
    }
}
