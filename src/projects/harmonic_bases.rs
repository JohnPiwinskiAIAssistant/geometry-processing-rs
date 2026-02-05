use crate::core::geometry::Geometry;
use crate::linear_algebra::DenseMatrix;
use crate::projects::vector_field_decomposition::HodgeDecomposition;

pub struct HarmonicBases<'a> {
    pub geometry: &'a Geometry<'a>,
}

impl<'a> HarmonicBases<'a> {
    pub fn new(geometry: &'a Geometry<'a>) -> Self {
        Self { geometry }
    }

    fn build_closed_primal_one_form(&self, generator: &[usize]) -> DenseMatrix {
        let e_count = self.geometry.mesh.edges.len();
        let mut omega = DenseMatrix::zeros(e_count, 1);
        for &h_idx in generator {
            let e_idx = self.geometry.mesh.halfedges[h_idx].edge.expect("Halfedge should have an edge");
            let sign = if self.geometry.mesh.edges[e_idx].halfedge == Some(h_idx) { 1.0 } else { -1.0 };
            omega.set(sign, e_idx, 0);
        }
        omega
    }

    pub fn compute(&self, hodge_decomposition: &HodgeDecomposition) -> Vec<DenseMatrix> {
        let mut gammas = Vec::new();
        for generator in &self.geometry.mesh.generators {
            let omega = self.build_closed_primal_one_form(generator);
            let d_alpha = hodge_decomposition.compute_exact_component(&omega);
            gammas.push(omega.minus(&d_alpha));
        }
        gammas
    }
}
