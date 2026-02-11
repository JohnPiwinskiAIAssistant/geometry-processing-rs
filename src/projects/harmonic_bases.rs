use crate::core::geometry::Geometry;
use crate::core::mesh::MeshBackend;
use crate::linear_algebra::DenseMatrix;
use crate::projects::vector_field_decomposition::HodgeDecomposition;

pub struct HarmonicBases<'a, B: MeshBackend> {
    pub geometry: &'a Geometry<'a, B>,
}

impl<'a, B: MeshBackend> HarmonicBases<'a, B> {
    pub fn new(geometry: &'a Geometry<'a, B>) -> Self {
        Self { geometry }
    }

    fn build_closed_primal_one_form(&self, generator: &[usize]) -> DenseMatrix {
        let e_count = self.geometry.mesh.num_edges();
        let mut omega = DenseMatrix::zeros(e_count, 1);
        for &h_idx in generator {
            let e_idx = self.geometry.mesh.backend.halfedge_edge(h_idx).expect("Halfedge should have an edge");
            let sign = if self.geometry.mesh.backend.edge_halfedge(e_idx) == Some(h_idx) { 1.0 } else { -1.0 };
            omega[(e_idx, 0)] = sign;
        }
        omega
    }

    pub fn compute(&self, hodge_decomposition: &HodgeDecomposition<B>) -> Vec<DenseMatrix> {
        let mut gammas = Vec::new();
        let num_generators = self.geometry.mesh.backend.num_generators();
        for i in 0..num_generators {
            let generator = self.geometry.mesh.backend.generator(i);
            let omega = self.build_closed_primal_one_form(generator);
            let d_alpha = hodge_decomposition.compute_exact_component(&omega);
            gammas.push(omega - d_alpha);
        }
        gammas
    }
}
