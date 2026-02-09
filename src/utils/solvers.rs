use crate::linear_algebra::{DenseMatrix, ComplexDenseMatrix, SparseMatrix, Complex};
use faer::prelude::SpSolver;
use num_complex::Complex64;

pub struct Solvers;

impl Solvers {
    pub fn residual(a: &SparseMatrix<Complex>, x: &ComplexDenseMatrix) -> f64 {
        let mut b = faer::Mat::<Complex64>::zeros(a.nrows(), x.ncols());
        let a_values = a.values();
        for j in 0..x.ncols() {
            for i_col in 0..a.ncols() {
                let start = a.col_ptrs()[i_col];
                let end = a.col_ptrs()[i_col+1];
                let xi = x.read(i_col, j);
                for k in start..end {
                    let mut val = b.read(a.row_indices()[k], j);
                    let a_val = Complex64::new(a_values.re[k], a_values.im[k]);
                    val += a_val * xi;
                    b.write(a.row_indices()[k], j, val);
                }
            }
        }
        b.norm_l1()
    }

    pub fn solve_inverse_power_method(a: &SparseMatrix<Complex>) -> ComplexDenseMatrix {
        let n = a.nrows();
        let mut x = faer::Mat::<Complex64>::from_fn(n, 1, |_, _| Complex64::new(rand::random(), rand::random()));
        
        let lu = a.as_ref().sp_lu().unwrap();
        
        for _ in 0..200 {
            let mut y_faer = x.clone();
            lu.solve_in_place(y_faer.as_mut());
            
            let mut norm_sq = 0.0;
            for i in 0..n {
                norm_sq += y_faer.read(i, 0).norm_sqr();
            }
            let norm = norm_sq.sqrt();
            
            for i in 0..n {
                x.write(i, 0, y_faer.read(i, 0) / norm);
            }
        }
        
        x
    }

    pub fn invert_2x2(m: &DenseMatrix) -> DenseMatrix {
        let a = m[(0, 0)];
        let b = m[(0, 1)];
        let c = m[(1, 0)];
        let d = m[(1, 1)];
        let det = a * d - b * c;
        let mut inv = DenseMatrix::zeros(2, 2);
        inv[(0, 0)] = d / det;
        inv[(0, 1)] = -b / det;
        inv[(1, 0)] = -c / det;
        inv[(1, 1)] = a / det;
        inv
    }
}
