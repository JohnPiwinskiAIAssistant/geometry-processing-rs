use crate::core::mesh::{Mesh, MeshBackend};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Vertex {
    pub index: usize,
}

impl Vertex {
    pub fn new(index: usize) -> Self {
        Self { index }
    }

    pub fn is_isolated<B: MeshBackend>(&self, mesh: &Mesh<B>) -> bool {
        mesh.is_isolated(self.index)
    }

    pub fn on_boundary<B: MeshBackend>(&self, mesh: &Mesh<B>) -> bool {
        mesh.on_boundary(self.index)
    }

    pub fn degree<B: MeshBackend>(&self, mesh: &Mesh<B>) -> usize {
        mesh.vertex_degree(self.index)
    }

    pub fn adjacent_vertices<'a, B: MeshBackend>(&self, mesh: &'a Mesh<B>, ccw: bool) -> Box<dyn Iterator<Item = usize> + 'a> {
        mesh.vertex_adjacent_vertices(self.index, ccw)
    }

    pub fn adjacent_edges<'a, B: MeshBackend>(&self, mesh: &'a Mesh<B>, ccw: bool) -> Box<dyn Iterator<Item = usize> + 'a> {
        mesh.vertex_adjacent_edges(self.index, ccw)
    }

    pub fn adjacent_faces<'a, B: MeshBackend>(&self, mesh: &'a Mesh<B>, ccw: bool) -> Box<dyn Iterator<Item = usize> + 'a> {
        mesh.vertex_adjacent_faces(self.index, ccw)
    }

    pub fn adjacent_halfedges<'a, B: MeshBackend>(&self, mesh: &'a Mesh<B>, ccw: bool) -> Box<dyn Iterator<Item = usize> + 'a> {
        mesh.vertex_adjacent_halfedges(self.index, ccw)
    }

    pub fn adjacent_corners<'a, B: MeshBackend>(&self, mesh: &'a Mesh<B>, ccw: bool) -> Box<dyn Iterator<Item = usize> + 'a> {
        mesh.vertex_adjacent_corners(self.index, ccw)
    }
}
