use crate::core::mesh::Mesh;
use crate::core::mesh_subset::MeshSubset;
use crate::linear_algebra::{SparseMatrix, Triplet, DenseMatrix};

pub struct SimplicialComplexOperators {
    pub mesh: Mesh,
    pub a0: SparseMatrix,
    pub a1: SparseMatrix,
}

impl SimplicialComplexOperators {
    pub fn new(mesh: Mesh) -> Self {
        let a0 = Self::build_vertex_edge_adjacency_matrix(&mesh);
        let a1 = Self::build_edge_face_adjacency_matrix(&mesh);
        Self { mesh, a0, a1 }
    }

    fn build_vertex_edge_adjacency_matrix(mesh: &Mesh) -> SparseMatrix {
        let mut t = Triplet::new(mesh.edges.len(), mesh.vertices.len());
        for e in &mesh.edges {
            let h_idx = e.halfedge.expect("Edge should have a halfedge");
            let v1 = mesh.halfedges[h_idx].vertex.expect("Halfedge should have a vertex");
            let next_idx = mesh.halfedges[h_idx].next.expect("Halfedge should have a next");
            let v2 = mesh.halfedges[next_idx].vertex.expect("Next should have a vertex");
            
            t.add_entry(1.0, e.index, v1);
            t.add_entry(1.0, e.index, v2);
        }
        SparseMatrix::from_triplet(t)
    }

    fn build_edge_face_adjacency_matrix(mesh: &Mesh) -> SparseMatrix {
        let mut t = Triplet::new(mesh.faces.len(), mesh.edges.len());
        for f in &mesh.faces {
            for h_idx in mesh.face_adjacent_halfedges(f.index, true) {
                let e_idx = mesh.halfedges[h_idx].edge.expect("Halfedge should have an edge");
                t.add_entry(1.0, f.index, e_idx);
            }
        }
        SparseMatrix::from_triplet(t)
    }

    pub fn build_vertex_vector(&self, subset: &MeshSubset) -> DenseMatrix {
        let mut v = DenseMatrix::zeros(self.mesh.vertices.len(), 1);
        for &v_idx in &subset.vertices {
            v.set(1.0, v_idx, 0);
        }
        v
    }

    pub fn build_edge_vector(&self, subset: &MeshSubset) -> DenseMatrix {
        let mut v = DenseMatrix::zeros(self.mesh.edges.len(), 1);
        for &e_idx in &subset.edges {
            v.set(1.0, e_idx, 0);
        }
        v
    }

    pub fn build_face_vector(&self, subset: &MeshSubset) -> DenseMatrix {
        let mut v = DenseMatrix::zeros(self.mesh.faces.len(), 1);
        for &f_idx in &subset.faces {
            v.set(1.0, f_idx, 0);
        }
        v
    }

    pub fn star(&self, subset: &MeshSubset) -> MeshSubset {
        let mut res = subset.clone();
        
        let edge_vector = self.a0.times_dense(&self.build_vertex_vector(subset));
        for i in 0..self.mesh.edges.len() {
            if edge_vector.get(i, 0) != 0.0 {
                res.add_edge(i);
            }
        }

        let face_vector = self.a1.times_dense(&self.build_edge_vector(&res));
        for i in 0..self.mesh.faces.len() {
            if face_vector.get(i, 0) != 0.0 {
                res.add_face(i);
            }
        }

        res
    }

    pub fn closure(&self, subset: &MeshSubset) -> MeshSubset {
        let mut res = subset.clone();
        
        let edge_vector = self.a1.transpose().times_dense(&self.build_face_vector(subset));
        for i in 0..self.mesh.edges.len() {
            if edge_vector.get(i, 0) != 0.0 {
                res.add_edge(i);
            }
        }

        let vertex_vector = self.a0.transpose().times_dense(&self.build_edge_vector(&res));
        for i in 0..self.mesh.vertices.len() {
            if vertex_vector.get(i, 0) != 0.0 {
                res.add_vertex(i);
            }
        }

        res
    }

    pub fn link(&self, subset: &MeshSubset) -> MeshSubset {
        let mut closure_star = self.closure(&self.star(subset));
        let star_closure = self.star(&self.closure(subset));

        closure_star.delete_subset(&star_closure);
        closure_star
    }

    pub fn is_complex(&self, subset: &MeshSubset) -> bool {
        subset.vertices == self.closure(subset).vertices &&
        subset.edges == self.closure(subset).edges &&
        subset.faces == self.closure(subset).faces
    }

    pub fn boundary(&self, subset: &MeshSubset) -> MeshSubset {
        let mut boundary = MeshSubset::new();
        if !subset.faces.is_empty() {
            let face_edges = self.a1.transpose().times_dense(&self.build_face_vector(subset));
            for i in 0..self.mesh.edges.len() {
                if face_edges.get(i, 0) == 1.0 {
                    boundary.add_edge(i);
                }
            }
        } else if !subset.edges.is_empty() {
            let edge_vertices = self.a0.transpose().times_dense(&self.build_edge_vector(subset));
            for i in 0..self.mesh.vertices.len() {
                if edge_vertices.get(i, 0) == 1.0 {
                    boundary.add_vertex(i);
                }
            }
        }
        self.closure(&boundary)
    }
}
