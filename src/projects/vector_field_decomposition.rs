use crate::core::geometry::Geometry;
use crate::core::dec::DEC;
use crate::linear_algebra::{SparseMatrix, DenseMatrix};

pub struct HodgeDecomposition {
    pub hodge1: SparseMatrix,
    pub hodge2: SparseMatrix,
    pub d0: SparseMatrix,
    pub d1: SparseMatrix,
    pub hodge1_inv: SparseMatrix,
    pub hodge2_inv: SparseMatrix,
    pub d0t: SparseMatrix,
    pub d1t: SparseMatrix,
    pub a: SparseMatrix,
    pub b: SparseMatrix,
}

impl HodgeDecomposition {
    pub fn new(geometry: &Geometry) -> Self {
        let hodge1 = DEC::build_hodge_star_1_form(geometry);
        let hodge2 = DEC::build_hodge_star_2_form(geometry);
        let d0 = DEC::build_exterior_derivative_0_form(geometry);
        let d1 = DEC::build_exterior_derivative_1_form(geometry);

        use crate::linear_algebra::sparse_matrix::SparseMatrixMethods;
        let hodge1_inv = hodge1.invert_diagonal();
        let hodge2_inv = hodge2.invert_diagonal();
        let d0t = d0.transpose().to_col_major().unwrap();
        let d1t = d1.transpose().to_col_major().unwrap();

        let v_count = geometry.mesh.vertices.len();
        let mut a = &d0t * &(&hodge1 * &d0);
        a = &a + &crate::linear_algebra::sparse_matrix::identity(v_count, v_count).scale(1e-8);

        let b = &d1 * &(&hodge1_inv * &d1t);

        Self {
            hodge1, hodge2, d0, d1,
            hodge1_inv, hodge2_inv,
            d0t, d1t,
            a, b,
        }
    }

    pub fn compute_exact_component(&self, omega: &DenseMatrix) -> DenseMatrix {
        let rhs = &self.d0t * &(&self.hodge1 * omega);
        let llt = crate::linear_algebra::Cholesky::new(&self.a);
        let alpha = llt.solve_positive_definite(&rhs);
        &self.d0 * &alpha
    }

    pub fn compute_co_exact_component(&self, omega: &DenseMatrix) -> DenseMatrix {
        let rhs = &self.d1 * omega;
        let lu = crate::linear_algebra::LU::new(&self.b);
        let beta_tilde = lu.solve_square(&rhs);
        &self.hodge1_inv * &(&self.d1t * &beta_tilde)
    }

    pub fn compute_harmonic_component(&self, omega: &DenseMatrix, d_alpha: &DenseMatrix, delta_beta: &DenseMatrix) -> DenseMatrix {
        omega - (d_alpha + delta_beta)
    }
}
