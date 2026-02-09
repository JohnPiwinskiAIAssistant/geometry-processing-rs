use faer::sparse::SparseColMat;
use crate::linear_algebra::{DenseMatrix};
use faer::prelude::SpSolver;
use num_complex::Complex64;

pub type SparseMatrix = SparseColMat<usize, f64>;
pub type ComplexSparseMatrix = SparseColMat<usize, Complex64>;

pub struct Triplet {
    pub rows: usize,
    pub cols: usize,
    pub data: Vec<(usize, usize, f64)>,
}

impl Triplet {
    pub fn new(rows: usize, cols: usize) -> Self {
        Self {
            rows,
            cols,
            data: Vec::new(),
        }
    }

    pub fn add_entry(&mut self, value: f64, row: usize, col: usize) {
        self.data.push((row, col, value));
    }
}

pub struct ComplexTriplet {
    pub rows: usize,
    pub cols: usize,
    pub data: Vec<(usize, usize, Complex64)>,
}

impl ComplexTriplet {
    pub fn new(rows: usize, cols: usize) -> Self {
        Self {
            rows,
            cols,
            data: Vec::new(),
        }
    }

    pub fn add_entry(&mut self, value: Complex64, row: usize, col: usize) {
        self.data.push((row, col, value));
    }
}

pub trait SparseMatrixMethods<E: faer::Entity> {
    fn from_triplets(rows: usize, cols: usize, triplets: &[(usize, usize, E)]) -> Self;
    fn invert_diagonal(&self) -> Self;
    fn sub_matrix(&self, r0: usize, r1: usize, c0: usize, c1: usize) -> Self;
    fn get_val(&self, row: usize, col: usize) -> E;
    fn frobenius_norm(&self) -> f64;
    fn scale(&self, s: E) -> Self;
    fn compute_nnz(&self) -> usize;
}

impl SparseMatrixMethods<f64> for SparseMatrix {
    fn from_triplets(rows: usize, cols: usize, triplets: &[(usize, usize, f64)]) -> Self {
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
                    data.push((i, j, 1.0 / self.values()[k]));
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
                    data.push((i - r0, j - c0, self.values()[k]));
                }
            }
        }
        Self::from_triplets(r1 - r0, c1 - c0, &data)
    }

    fn get_val(&self, row: usize, col: usize) -> f64 {
        if col >= self.ncols() { return 0.0; }
        let start = self.col_ptrs()[col];
        let end = self.col_ptrs()[col+1];
        let row_indices = self.row_indices();
        let values = self.values();
        for k in start..end {
            if row_indices[k] == row {
                return values[k];
            }
        }
        0.0
    }

    fn frobenius_norm(&self) -> f64 {
        let mut sum_sq = 0.0;
        for &val in self.values() {
            sum_sq += val * val;
        }
        sum_sq.sqrt()
    }

    fn scale(&self, s: f64) -> Self {
        let mut data = Vec::new();
        for j in 0..self.ncols() {
            let start = self.col_ptrs()[j];
            let end = self.col_ptrs()[j+1];
            for k in start..end {
                data.push((self.row_indices()[k], j, self.values()[k] * s));
            }
        }
        Self::from_triplets(self.nrows(), self.ncols(), &data)
    }

    fn compute_nnz(&self) -> usize {
        self.values().len()
    }
}

impl SparseMatrixMethods<Complex64> for ComplexSparseMatrix {
    fn from_triplets(rows: usize, cols: usize, triplets: &[(usize, usize, Complex64)]) -> Self {
        SparseColMat::try_new_from_triplets(rows, cols, triplets).unwrap()
    }

    fn invert_diagonal(&self) -> Self {
        let mut data = Vec::new();
        let values = self.values();
        for j in 0..self.ncols() {
            let start = self.col_ptrs()[j];
            let end = self.col_ptrs()[j+1];
            for k in start..end {
                let i = self.row_indices()[k];
                if i == j {
                    let val = Complex64::new(values.re[k], values.im[k]);
                    data.push((i, j, 1.0 / val));
                }
            }
        }
        Self::from_triplets(self.nrows(), self.ncols(), &data)
    }

    fn sub_matrix(&self, r0: usize, r1: usize, c0: usize, c1: usize) -> Self {
        let mut data = Vec::new();
        let values = self.values();
        for j in c0..c1 {
            let start = self.col_ptrs()[j];
            let end = self.col_ptrs()[j+1];
            for k in start..end {
                let i = self.row_indices()[k];
                if i >= r0 && i < r1 {
                    let val = Complex64::new(values.re[k], values.im[k]);
                    data.push((i - r0, j - c0, val));
                }
            }
        }
        Self::from_triplets(r1 - r0, c1 - c0, &data)
    }

    fn get_val(&self, row: usize, col: usize) -> Complex64 {
        if col >= self.ncols() { return Complex64::new(0.0, 0.0); }
        let start = self.col_ptrs()[col];
        let end = self.col_ptrs()[col+1];
        let row_indices = self.row_indices();
        let values = self.values();
        for k in start..end {
            if row_indices[k] == row {
                return Complex64::new(values.re[k], values.im[k]);
            }
        }
        Complex64::new(0.0, 0.0)
    }

    fn frobenius_norm(&self) -> f64 {
        let mut sum_sq = 0.0;
        let values = self.values();
        for i in 0..values.re.len() {
            sum_sq += values.re[i] * values.re[i] + values.im[i] * values.im[i];
        }
        sum_sq.sqrt()
    }

    fn scale(&self, s: Complex64) -> Self {
        let mut data = Vec::new();
        let values = self.values();
        for j in 0..self.ncols() {
            let start = self.col_ptrs()[j];
            let end = self.col_ptrs()[j+1];
            for k in start..end {
                let val = Complex64::new(values.re[k], values.im[k]);
                data.push((self.row_indices()[k], j, val * s));
            }
        }
        Self::from_triplets(self.nrows(), self.ncols(), &data)
    }

    fn compute_nnz(&self) -> usize {
        self.values().re.len()
    }
}

pub fn diag(v: &DenseMatrix) -> SparseMatrix {
    let n = v.nrows();
    let data: Vec<_> = (0..n).map(|i| (i, i, v[(i, 0)])).collect();
    SparseMatrix::from_triplets(n, n, &data)
}

pub fn identity(m: usize, n: usize) -> SparseMatrix {
    let data: Vec<_> = (0..m.min(n)).map(|i| (i, i, 1.0)).collect();
    SparseMatrix::from_triplets(m, n, &data)
}

pub struct Cholesky {
    pub llt: faer::sparse::linalg::solvers::Cholesky<usize, f64>,
}

impl Cholesky {
    pub fn new(mat: &SparseMatrix) -> Self {
        Self { llt: mat.as_ref().sp_cholesky(faer::Side::Lower).unwrap() }
    }
    pub fn solve_positive_definite(&self, rhs: &DenseMatrix) -> DenseMatrix {
        let mut res = rhs.clone();
        self.llt.solve_in_place(res.as_mut());
        res
    }
}

pub struct LU {
    pub lu: faer::sparse::linalg::solvers::Lu<usize, f64>,
}

impl LU {
    pub fn new(mat: &SparseMatrix) -> Self {
        Self { lu: mat.as_ref().sp_lu().unwrap() }
    }
    pub fn solve_square(&self, rhs: &DenseMatrix) -> DenseMatrix {
        let mut res = rhs.clone();
        self.lu.solve_in_place(res.as_mut());
        res
    }
}

pub struct QR {
    pub mat: SparseMatrix,
}

impl QR {
    pub fn new(mat: &SparseMatrix) -> Self {
        Self { mat: mat.clone() }
    }
    pub fn solve(&self, rhs: &DenseMatrix) -> DenseMatrix {
        let qr = self.mat.to_dense().as_ref().col_piv_qr();
        let mut res = rhs.clone();
        qr.solve_in_place(res.as_mut());
        res
    }
}
