use faer::sparse::SparseColMat;
use faer::prelude::SpSolver;
use num_complex::Complex64;
use crate::linear_algebra::traits::{Scalar, SparseOps, LinearSolver};
use crate::linear_algebra::{DenseMatrix};

pub type SparseMatrix<S> = SparseColMat<usize, S>;

pub struct Triplet<S: Scalar> {
    pub rows: usize,
    pub cols: usize,
    pub data: Vec<(usize, usize, S)>,
}

impl<S: Scalar> Triplet<S> {
    pub fn new(rows: usize, cols: usize) -> Self {
        Self {
            rows,
            cols,
            data: Vec::new(),
        }
    }

    pub fn add_entry(&mut self, value: S, row: usize, col: usize) {
        self.data.push((row, col, value));
    }
}

impl<S: Scalar> SparseOps<S> for SparseMatrix<S> {
    fn from_triplets(rows: usize, cols: usize, triplets: &[(usize, usize, S)]) -> Self {
        SparseColMat::try_new_from_triplets(rows, cols, triplets).unwrap()
    }

    fn invert_diagonal(&self) -> Self {
        let mut data = Vec::new();
        for j in 0..self.ncols() {
            let start = self.col_ptrs()[j];
            let end = self.col_ptrs()[j+1];
            for k in start..end {
                let i = self.row_indices()[k];
                if i == j {
                    data.push((i, j, S::one() / self.get_val(i, j)));
                }
            }
        }
        Self::from_triplets(self.nrows(), self.ncols(), &data)
    }

    fn sub_matrix(&self, r0: usize, r1: usize, c0: usize, c1: usize) -> Self {
        let mut data = Vec::new();
        for j in c0..c1 {
            let start = self.col_ptrs()[j];
            let end = self.col_ptrs()[j+1];
            for k in start..end {
                let i = self.row_indices()[k];
                if i >= r0 && i < r1 {
                    data.push((i - r0, j - c0, self.get_val(i, j)));
                }
            }
        }
        Self::from_triplets(r1 - r0, c1 - c0, &data)
    }

    fn get_val(&self, row: usize, col: usize) -> S {
        if col >= self.ncols() { return S::zero(); }
        let start = self.col_ptrs()[col];
        let end = self.col_ptrs()[col+1];
        let row_indices = self.row_indices();
        
        for k in start..end {
            if row_indices[k] == row {
                return S::get_sparse_value(self, k);
            }
        }
        S::zero()
    }

    fn frobenius_norm(&self) -> f64 {
        let mut sum_sq = 0.0;
        for j in 0..self.ncols() {
            let start = self.col_ptrs()[j];
            let end = self.col_ptrs()[j+1];
            for k in start..end {
                let val = S::get_sparse_value(self, k);
                sum_sq += val.abs() * val.abs();
            }
        }
        sum_sq.sqrt()
    }

    fn scale(&self, s: S) -> Self {
        let mut data = Vec::new();
        for j in 0..self.ncols() {
            let start = self.col_ptrs()[j];
            let end = self.col_ptrs()[j+1];
            for k in start..end {
                let val = S::get_sparse_value(self, k);
                data.push((self.row_indices()[k], j, val * s));
            }
        }
        Self::from_triplets(self.nrows(), self.ncols(), &data)
    }

    fn compute_nnz(&self) -> usize {
        S::sparse_values_count(self)
    }
}

pub fn diag<S: Scalar>(v: &faer::Mat<S>) -> SparseMatrix<S> {
    let n = v.nrows();
    let mut triplet = Triplet::new(n, n);
    for i in 0..n {
        triplet.add_entry(v.read(i, 0), i, i);
    }
    SparseOps::from_triplets(n, n, &triplet.data)
}

pub fn identity<S: Scalar>(m: usize, n: usize) -> SparseMatrix<S> {
    let mut triplet = Triplet::new(m, n);
    for i in 0..m.min(n) {
        triplet.add_entry(S::one(), i, i);
    }
    SparseOps::from_triplets(m, n, &triplet.data)
}

pub struct Cholesky<S: Scalar> {
    pub llt: faer::sparse::linalg::solvers::Cholesky<usize, S>,
}

impl<S: Scalar> Cholesky<S> {
    pub fn new(mat: &SparseMatrix<S>) -> Self {
        Self { llt: mat.as_ref().sp_cholesky(faer::Side::Lower).unwrap() }
    }
}

impl<S: Scalar> LinearSolver<S> for Cholesky<S> {
    fn solve(&self, rhs: &faer::Mat<S>) -> faer::Mat<S> {
        let mut res = rhs.clone();
        self.llt.solve_in_place(res.as_mut());
        res
    }
}

pub struct LU<S: Scalar> {
    pub lu: faer::sparse::linalg::solvers::Lu<usize, S>,
}

impl<S: Scalar> LU<S> {
    pub fn new(mat: &SparseMatrix<S>) -> Self {
        Self { lu: mat.as_ref().sp_lu().unwrap() }
    }
}

impl<S: Scalar> LinearSolver<S> for LU<S> {
    fn solve(&self, rhs: &faer::Mat<S>) -> faer::Mat<S> {
        let mut res = rhs.clone();
        self.lu.solve_in_place(res.as_mut());
        res
    }
}

pub struct QR<S: Scalar> {
    pub mat: SparseMatrix<S>,
}

impl<S: Scalar> QR<S> {
    pub fn new(mat: &SparseMatrix<S>) -> Self {
        Self { mat: mat.clone() }
    }
}

impl<S: Scalar> LinearSolver<S> for QR<S> {
    fn solve(&self, rhs: &faer::Mat<S>) -> faer::Mat<S> {
        // QR is currently not implemented in faer for sparse, we use dense fallback
        let qr = self.mat.to_dense().as_ref().col_piv_qr();
        let mut res = rhs.clone();
        qr.solve_in_place(res.as_mut());
        res
    }
}
