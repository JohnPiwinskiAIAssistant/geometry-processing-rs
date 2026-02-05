use std::collections::HashSet;

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct MeshSubset {
    pub vertices: HashSet<usize>,
    pub edges: HashSet<usize>,
    pub faces: HashSet<usize>,
}

impl MeshSubset {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn with_elements(vertices: HashSet<usize>, edges: HashSet<usize>, faces: HashSet<usize>) -> Self {
        Self { vertices, edges, faces }
    }

    pub fn deep_copy(&self) -> Self {
        self.clone()
    }

    pub fn reset(&mut self) {
        self.vertices.clear();
        self.edges.clear();
        self.faces.clear();
    }

    pub fn add_vertex(&mut self, v: usize) { self.vertices.insert(v); }
    pub fn delete_vertex(&mut self, v: usize) { self.vertices.remove(&v); }

    pub fn add_edge(&mut self, e: usize) { self.edges.insert(e); }
    pub fn delete_edge(&mut self, e: usize) { self.edges.remove(&e); }

    pub fn add_face(&mut self, f: usize) { self.faces.insert(f); }
    pub fn delete_face(&mut self, f: usize) { self.faces.remove(&f); }

    pub fn add_subset(&mut self, other: &MeshSubset) {
        for &v in &other.vertices { self.add_vertex(v); }
        for &e in &other.edges { self.add_edge(e); }
        for &f in &other.faces { self.add_face(f); }
    }

    pub fn delete_subset(&mut self, other: &MeshSubset) {
        for &v in &other.vertices { self.delete_vertex(v); }
        for &e in &other.edges { self.delete_edge(e); }
        for &f in &other.faces { self.delete_face(f); }
    }

    pub fn equals(&self, other: &MeshSubset) -> bool {
        self == other
    }
}
