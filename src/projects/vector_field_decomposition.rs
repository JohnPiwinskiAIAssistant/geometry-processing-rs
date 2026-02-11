use crate::core::geometry::Geometry;
use crate::core::mesh::MeshBackend;
use crate::core::dec::DEC;
use crate::linear_algebra::{SparseMatrix, DenseMatrix, Cholesky, LU};
use crate::linear_algebra::traits::{SparseOps, LinearSolver};

pub struct HodgeDecomposition<B: MeshBackend> {
    pub hodge1: SparseMatrix<f64>,
    pub hodge2: SparseMatrix<f64>,
    pub d0: SparseMatrix<f64>,
    pub d1: SparseMatrix<f64>,
    pub hodge1_inv: SparseMatrix<f64>,
    pub hodge2_inv: SparseMatrix<f64>,
    pub d0t: SparseMatrix<f64>,
    pub d1t: SparseMatrix<f64>,
    pub a: SparseMatrix<f64>,
    pub b: SparseMatrix<f64>,
    pub _marker: std::marker::PhantomData<B>,
}

impl<B: MeshBackend> HodgeDecomposition<B> {
    pub fn new(geometry: &Geometry<B>) -> Self {
        let hodge1 = DEC::build_hodge_star_1_form(geometry);
        let hodge2 = DEC::build_hodge_star_2_form(geometry);
        let d0 = DEC::build_exterior_derivative_0_form(geometry);
        let d1 = DEC::build_exterior_derivative_1_form(geometry);

        let hodge1_inv = hodge1.invert_diagonal();
        let hodge2_inv = hodge2.invert_diagonal();
        let d0t = d0.transpose().to_col_major().unwrap();
        let d1t = d1.transpose().to_col_major().unwrap();

        let v_count = geometry.mesh.num_vertices();
        let mut a = &d0t * &(&hodge1 * &d0);
        a = &a + &crate::linear_algebra::sparse_matrix::identity::<f64>(v_count, v_count).scale(1e-8);

        let b = &d1 * &(&hodge1_inv * &d1t);

        Self {
            hodge1, hodge2, d0, d1,
            hodge1_inv, hodge2_inv,
            d0t, d1t,
            a, b,
            _marker: std::marker::PhantomData,
        }
    }

    pub fn compute_exact_component(&self, omega: &DenseMatrix) -> DenseMatrix {
        let rhs = &self.d0t * &(&self.hodge1 * omega);
        let llt = Cholesky::<f64>::new(&self.a);
        let alpha = llt.solve(&rhs);
        &self.d0 * &alpha
    }

    pub fn compute_co_exact_component(&self, omega: &DenseMatrix) -> DenseMatrix {
        let rhs = &self.d1 * omega;
        let lu = LU::<f64>::new(&self.b);
        let beta_tilde = lu.solve(&rhs);
        &self.hodge1_inv * &(&self.d1t * &beta_tilde)
    }

    pub fn compute_harmonic_component(&self, omega: &DenseMatrix, d_alpha: &DenseMatrix, delta_beta: &DenseMatrix) -> DenseMatrix {
        omega - (d_alpha + delta_beta)
    }
}
