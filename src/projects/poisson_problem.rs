use crate::core::geometry::Geometry;
use crate::linear_algebra::{DenseMatrix, SparseMatrix};

pub struct ScalarPoissonProblem<'a> {
    pub geometry: &'a Geometry<'a>,
    pub a: SparseMatrix,
    pub m: SparseMatrix,
    pub total_area: f64,
}

impl<'a> ScalarPoissonProblem<'a> {
    pub fn new(geometry: &'a Geometry<'a>) -> Self {
        let a = geometry.laplace_matrix();
        let m = geometry.mass_matrix();
        let total_area = geometry.total_area();
        Self { geometry, a, m, total_area }
    }

    pub fn solve(&self, rho: &DenseMatrix) -> DenseMatrix {
        let v_count = self.m.n_rows();
        let total_rho = self.m.times_dense(rho).sum();
        
        let rho_bar_val = total_rho / self.total_area;
        let rho_bar = DenseMatrix::constant(rho_bar_val, v_count, 1);
        
        // rhs = M * (rho_bar - rho)
        let rhs = self.m.times_dense(&rho_bar.minus(rho));

        let llt = self.a.chol();
        llt.solve_positive_definite(&rhs)
    }
}
