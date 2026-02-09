use crate::core::geometry::Geometry;
use crate::linear_algebra::{DenseMatrix, SparseMatrix, Cholesky};
use crate::linear_algebra::traits::LinearSolver;

pub struct ScalarPoissonProblem<'a> {
    pub geometry: &'a Geometry<'a>,
    pub a: SparseMatrix<f64>,
    pub m: SparseMatrix<f64>,
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
        let v_count = self.m.nrows();
        let total_rho = (&self.m * rho).sum();
        
        let rho_bar_val = total_rho / self.total_area;
        let rho_bar = DenseMatrix::from_fn(v_count, 1, |_, _| rho_bar_val);
        
        // rhs = M * (rho_bar - rho)
        let rhs = &self.m * &(rho_bar - rho);

        let llt = Cholesky::<f64>::new(&self.a);
        llt.solve(&rhs)
    }
}
