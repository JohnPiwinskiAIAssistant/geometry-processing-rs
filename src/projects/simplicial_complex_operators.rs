use crate::core::mesh::{Mesh, MeshBackend};
use crate::core::mesh_subset::MeshSubset;
use crate::linear_algebra::{SparseMatrix, Triplet, DenseMatrix};
use crate::linear_algebra::traits::SparseOps;

pub struct SimplicialComplexOperators<B: MeshBackend> {
    pub mesh: Mesh<B>,
    pub a0: SparseMatrix<f64>,
    pub a1: SparseMatrix<f64>,
}

impl<B: MeshBackend> SimplicialComplexOperators<B> {
    pub fn new(mesh: Mesh<B>) -> Self {
        let a0 = Self::build_vertex_edge_adjacency_matrix(&mesh);
        let a1 = Self::build_edge_face_adjacency_matrix(&mesh);
        Self { mesh, a0, a1 }
    }

    fn build_vertex_edge_adjacency_matrix(mesh: &Mesh<B>) -> SparseMatrix<f64> {
        let num_edges = mesh.num_edges();
        let num_vertices = mesh.num_vertices();
        let mut t = Triplet::<f64>::new(num_edges, num_vertices);
        for i in 0..num_edges {
            let h_idx = mesh.backend.edge_halfedge(i).expect("Edge should have a halfedge");
            let v1 = mesh.backend.halfedge_vertex(h_idx).expect("Halfedge should have a vertex");
            let next_idx = mesh.backend.halfedge_next(h_idx).expect("Halfedge should have a next");
            let v2 = mesh.backend.halfedge_vertex(next_idx).expect("Next should have a vertex");
            
            t.add_entry(1.0, i, v1);
            t.add_entry(1.0, i, v2);
        }
        SparseOps::from_triplets(num_edges, num_vertices, &t.data)
    }

    fn build_edge_face_adjacency_matrix(mesh: &Mesh<B>) -> SparseMatrix<f64> {
        let num_faces = mesh.num_faces();
        let num_edges = mesh.num_edges();
        let mut t = Triplet::<f64>::new(num_faces, num_edges);
        for f_idx in 0..num_faces {
            for h_idx in mesh.face_adjacent_halfedges(f_idx, true) {
                let e_idx = mesh.backend.halfedge_edge(h_idx).expect("Halfedge should have an edge");
                t.add_entry(1.0, f_idx, e_idx);
            }
        }
        SparseOps::from_triplets(num_faces, num_edges, &t.data)
    }

    pub fn build_vertex_vector(&self, subset: &MeshSubset) -> DenseMatrix {
        let mut v = DenseMatrix::zeros(self.mesh.num_vertices(), 1);
        for &v_idx in &subset.vertices {
            v[(v_idx, 0)] = 1.0;
        }
        v
    }

    pub fn build_edge_vector(&self, subset: &MeshSubset) -> DenseMatrix {
        let mut v = DenseMatrix::zeros(self.mesh.num_edges(), 1);
        for &e_idx in &subset.edges {
            v[(e_idx, 0)] = 1.0;
        }
        v
    }

    pub fn build_face_vector(&self, subset: &MeshSubset) -> DenseMatrix {
        let mut v = DenseMatrix::zeros(self.mesh.num_faces(), 1);
        for &f_idx in &subset.faces {
            v[(f_idx, 0)] = 1.0;
        }
        v
    }

    pub fn star(&self, subset: &MeshSubset) -> MeshSubset {
        let mut res = subset.clone();
        
        let edge_vector = &self.a0 * &self.build_vertex_vector(subset);
        let num_edges = self.mesh.num_edges();
        for i in 0..num_edges {
            if edge_vector[(i, 0)] != 0.0 {
                res.add_edge(i);
            }
        }

        let face_vector = &self.a1 * &self.build_edge_vector(&res);
        let num_faces = self.mesh.num_faces();
        for i in 0..num_faces {
            if face_vector[(i, 0)] != 0.0 {
                res.add_face(i);
            }
        }

        res
    }

    pub fn closure(&self, subset: &MeshSubset) -> MeshSubset {
        let mut res = subset.clone();
        
        let edge_vector = &self.a1.transpose().to_col_major().unwrap() * &self.build_face_vector(subset);
        let num_edges = self.mesh.num_edges();
        for i in 0..num_edges {
            if edge_vector[(i, 0)] != 0.0 {
                res.add_edge(i);
            }
        }

        let vertex_vector = &self.a0.transpose().to_col_major().unwrap() * &self.build_edge_vector(&res);
        let num_vertices = self.mesh.num_vertices();
        for i in 0..num_vertices {
            if vertex_vector[(i, 0)] != 0.0 {
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
            let face_edges = &self.a1.transpose().to_col_major().unwrap() * &self.build_face_vector(subset);
            let num_edges = self.mesh.num_edges();
            for i in 0..num_edges {
                if face_edges[(i, 0)] == 1.0 {
                    boundary.add_edge(i);
                }
            }
        } else if !subset.edges.is_empty() {
            let edge_vertices = &self.a0.transpose().to_col_major().unwrap() * &self.build_edge_vector(subset);
            let num_vertices = self.mesh.num_vertices();
            for i in 0..num_vertices {
                if edge_vertices[(i, 0)] == 1.0 {
                    boundary.add_vertex(i);
                }
            }
        }
        self.closure(&boundary)
    }
}
