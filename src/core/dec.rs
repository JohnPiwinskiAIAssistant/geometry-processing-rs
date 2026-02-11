use crate::core::geometry::Geometry;
use crate::core::mesh::{MeshBackend, Vertex, Face};
use crate::linear_algebra::{SparseMatrix, Triplet};
use crate::linear_algebra::traits::SparseOps;

pub struct DEC;

impl DEC {
    pub fn build_hodge_star_0_form<B: MeshBackend>(geometry: &Geometry<B>) -> SparseMatrix<f64> {
        let v_count = geometry.mesh.num_vertices();
        let mut triplet = Triplet::<f64>::new(v_count, v_count);
        for i in 0..v_count {
            let v = Vertex::new(i);
            let area = geometry.barycentric_dual_area(&v);
            triplet.add_entry(area, i, i);
        }
        SparseOps::from_triplets(v_count, v_count, &triplet.data)
    }

    pub fn build_hodge_star_1_form<B: MeshBackend>(geometry: &Geometry<B>) -> SparseMatrix<f64> {
        let e_count = geometry.mesh.num_edges();
        let mut triplet = Triplet::<f64>::new(e_count, e_count);
        for i in 0..e_count {
            let h_idx = geometry.mesh.backend.edge_halfedge(i).expect("Edge should have a halfedge");
            let twin_idx = geometry.mesh.backend.halfedge_twin(h_idx).expect("Halfedge should have a twin");
            let w = (geometry.cotan(h_idx) + geometry.cotan(twin_idx)) / 2.0;
            triplet.add_entry(w, i, i);
        }
        SparseOps::from_triplets(e_count, e_count, &triplet.data)
    }

    pub fn build_hodge_star_2_form<B: MeshBackend>(geometry: &Geometry<B>) -> SparseMatrix<f64> {
        let f_count = geometry.mesh.num_faces();
        let mut triplet = Triplet::<f64>::new(f_count, f_count);
        for i in 0..f_count {
            let area = geometry.area(&Face::new(i));
            triplet.add_entry(1.0 / area, i, i);
        }
        SparseOps::from_triplets(f_count, f_count, &triplet.data)
    }

    pub fn build_exterior_derivative_0_form<B: MeshBackend>(geometry: &Geometry<B>) -> SparseMatrix<f64> {
        let e_count = geometry.mesh.num_edges();
        let v_count = geometry.mesh.num_vertices();
        let mut triplet = Triplet::<f64>::new(e_count, v_count);
        for i in 0..e_count {
            let h_idx = geometry.mesh.backend.edge_halfedge(i).expect("Edge should have a halfedge");
            let twin_idx = geometry.mesh.backend.halfedge_twin(h_idx).expect("Halfedge should have a twin");
            
            let j = geometry.mesh.backend.halfedge_vertex(h_idx).expect("Halfedge should have a vertex");
            let k = geometry.mesh.backend.halfedge_vertex(twin_idx).expect("Twin should have a vertex");

            triplet.add_entry(1.0, i, j);
            triplet.add_entry(-1.0, i, k);
        }
        SparseOps::from_triplets(e_count, v_count, &triplet.data)
    }

    pub fn build_exterior_derivative_1_form<B: MeshBackend>(geometry: &Geometry<B>) -> SparseMatrix<f64> {
        let f_count = geometry.mesh.num_faces();
        let e_count = geometry.mesh.num_edges();
        let mut triplet = Triplet::<f64>::new(f_count, e_count);
        for i in 0..f_count {
            for h_idx in geometry.mesh.face_adjacent_halfedges(i, true) {
                let e_idx = geometry.mesh.backend.halfedge_edge(h_idx).expect("Halfedge should have an edge");
                let sign = if geometry.mesh.backend.edge_halfedge(e_idx) == Some(h_idx) { 1.0 } else { -1.0 };
                triplet.add_entry(sign, i, e_idx);
            }
        }
        SparseOps::from_triplets(f_count, e_count, &triplet.data)
    }
}
