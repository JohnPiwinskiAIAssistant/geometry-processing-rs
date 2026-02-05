use crate::linear_algebra::{DenseMatrix, ComplexDenseMatrix, ComplexSparseMatrix, Complex};
use faer::prelude::SpSolverCore;
use num_complex::Complex64;

pub struct Solvers;

impl Solvers {
    pub fn residual(a: &ComplexSparseMatrix, x: &ComplexDenseMatrix) -> f64 {
        let mut b = ComplexDenseMatrix::zeros(a.n_rows(), x.n_cols());
        let a_values = a.mat.values();
        for j in 0..x.n_cols() {
            for i_col in 0..a.n_cols() {
                let start = a.mat.col_ptrs()[i_col];
                let end = a.mat.col_ptrs()[i_col+1];
                let xi = x.get(i_col, j);
                for k in start..end {
                    let mut val = b.get(a.mat.row_indices()[k], j);
                    let a_val = Complex64::new(a_values.re[k], a_values.im[k]);
                    val.inner += a_val * xi.inner;
                    b.set(val, a.mat.row_indices()[k], j);
                }
            }
        }
        b.norm(2)
    }

    pub fn solve_inverse_power_method(a: &ComplexSparseMatrix) -> ComplexDenseMatrix {
        let n = a.n_rows();
        let mut x = ComplexDenseMatrix::random(n, 1);
        
        let lu = a.mat.as_ref().sp_lu().unwrap();
        
        for _ in 0..200 {
            let mut y_faer = faer::Mat::<Complex64>::zeros(n, 1);
            for i in 0..n {
                y_faer.write(i, 0, x.get(i, 0).inner);
            }
            
            lu.solve_in_place_with_conj_impl(y_faer.as_mut(), faer::Conj::No);
            
            let mut norm_sq = 0.0;
            for i in 0..n {
                norm_sq += y_faer.read(i, 0).norm_sqr();
            }
            let norm = norm_sq.sqrt();
            
            for i in 0..n {
                x.set(Complex { inner: y_faer.read(i, 0) / norm }, i, 0);
            }
        }
        
        x
    }

    pub fn invert_2x2(m: &DenseMatrix) -> DenseMatrix {
        let a = m.get(0, 0);
        let b = m.get(0, 1);
        let c = m.get(1, 0);
        let d = m.get(1, 1);
        let det = a * d - b * c;
        let mut inv = DenseMatrix::zeros(2, 2);
        inv.set(d / det, 0, 0);
        inv.set(-b / det, 0, 1);
        inv.set(-c / det, 1, 0);
        inv.set(a / det, 1, 1);
        inv
    }
}
