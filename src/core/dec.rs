use crate::core::geometry::Geometry;
use crate::linear_algebra::{SparseMatrix, Triplet};
use crate::linear_algebra::traits::SparseOps;

pub struct DEC;

impl DEC {
    pub fn build_hodge_star_0_form(geometry: &Geometry) -> SparseMatrix<f64> {
        let v_count = geometry.mesh.vertices.len();
        let mut triplet = Triplet::<f64>::new(v_count, v_count);
        for v in &geometry.mesh.vertices {
            let area = geometry.barycentric_dual_area(v);
            triplet.add_entry(area, v.index, v.index);
        }
        SparseOps::from_triplets(v_count, v_count, &triplet.data)
    }

    pub fn build_hodge_star_1_form(geometry: &Geometry) -> SparseMatrix<f64> {
        let e_count = geometry.mesh.edges.len();
        let mut triplet = Triplet::<f64>::new(e_count, e_count);
        for e in &geometry.mesh.edges {
            let h_idx = e.halfedge.expect("Edge should have a halfedge");
            let twin_idx = geometry.mesh.halfedges[h_idx].twin.expect("Halfedge should have a twin");
            let w = (geometry.cotan(h_idx) + geometry.cotan(twin_idx)) / 2.0;
            triplet.add_entry(w, e.index, e.index);
        }
        SparseOps::from_triplets(e_count, e_count, &triplet.data)
    }

    pub fn build_hodge_star_2_form(geometry: &Geometry) -> SparseMatrix<f64> {
        let f_count = geometry.mesh.faces.len();
        let mut triplet = Triplet::<f64>::new(f_count, f_count);
        for f in &geometry.mesh.faces {
            let area = geometry.area(f);
            triplet.add_entry(1.0 / area, f.index, f.index);
        }
        SparseOps::from_triplets(f_count, f_count, &triplet.data)
    }

    pub fn build_exterior_derivative_0_form(geometry: &Geometry) -> SparseMatrix<f64> {
        let e_count = geometry.mesh.edges.len();
        let v_count = geometry.mesh.vertices.len();
        let mut triplet = Triplet::<f64>::new(e_count, v_count);
        for e in &geometry.mesh.edges {
            let h_idx = e.halfedge.expect("Edge should have a halfedge");
            let twin_idx = geometry.mesh.halfedges[h_idx].twin.expect("Halfedge should have a twin");
            
            let j = geometry.mesh.halfedges[h_idx].vertex.expect("Halfedge should have a vertex");
            let k = geometry.mesh.halfedges[twin_idx].vertex.expect("Twin should have a vertex");

            triplet.add_entry(1.0, e.index, j);
            triplet.add_entry(-1.0, e.index, k);
        }
        SparseOps::from_triplets(e_count, v_count, &triplet.data)
    }

    pub fn build_exterior_derivative_1_form(geometry: &Geometry) -> SparseMatrix<f64> {
        let f_count = geometry.mesh.faces.len();
        let e_count = geometry.mesh.edges.len();
        let mut triplet = Triplet::<f64>::new(f_count, e_count);
        for f in &geometry.mesh.faces {
            for h_idx in geometry.mesh.face_adjacent_halfedges(f.index, true) {
                let e_idx = geometry.mesh.halfedges[h_idx].edge.expect("Halfedge should have an edge");
                let sign = if geometry.mesh.edges[e_idx].halfedge == Some(h_idx) { 1.0 } else { -1.0 };
                triplet.add_entry(sign, f.index, e_idx);
            }
        }
        SparseOps::from_triplets(f_count, e_count, &triplet.data)
    }
}
