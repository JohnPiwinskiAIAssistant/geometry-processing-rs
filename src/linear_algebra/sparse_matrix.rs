use faer::sparse::SparseColMat;
use crate::linear_algebra::{Complex, DenseMatrix, ComplexDenseMatrix};
use num_complex::Complex64;
use faer::Mat;
use faer::prelude::SpSolverCore;

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

    pub fn add_entry(&mut self, value: Complex, row: usize, col: usize) {
        self.data.push((row, col, value.inner));
    }
}

pub struct SparseMatrix {
    pub mat: SparseColMat<usize, f64>,
}

impl SparseMatrix {
    pub fn from_triplet(triplet: Triplet) -> Self {
        let mat = SparseColMat::try_new_from_triplets(
            triplet.rows,
            triplet.cols,
            &triplet.data,
        ).unwrap();
        
        Self { mat }
    }

    pub fn get(&self, row: usize, col: usize) -> f64 {
        if col >= self.mat.ncols() { return 0.0; }
        let start = self.mat.col_ptrs()[col];
        let end = self.mat.col_ptrs()[col+1];
        let row_indices = self.mat.row_indices();
        let values = self.mat.values();
        for k in start..end {
            if row_indices[k] == row {
                return values[k];
            }
        }
        0.0
    }

    pub fn n_rows(&self) -> usize { self.mat.nrows() }
    pub fn n_cols(&self) -> usize { self.mat.ncols() }
    pub fn nnz(&self) -> usize { self.mat.compute_nnz() }

    pub fn frobenius_norm(&self) -> f64 {
        let mut sum_sq = 0.0;
        for &val in self.mat.values() {
            sum_sq += val * val;
        }
        sum_sq.sqrt()
    }

    pub fn to_dense(&self) -> DenseMatrix {
        DenseMatrix { inner: self.mat.to_dense() }
    }

    pub fn plus(&self, other: &SparseMatrix) -> SparseMatrix {
        let mut data = Vec::new();
        for i in 0..self.mat.ncols() {
            let start = self.mat.col_ptrs()[i];
            let end = self.mat.col_ptrs()[i+1];
            for k in start..end {
                data.push((self.mat.row_indices()[k], i, self.mat.values()[k]));
            }
        }
        for i in 0..other.mat.ncols() {
            let start = other.mat.col_ptrs()[i];
            let end = other.mat.col_ptrs()[i+1];
            for k in start..end {
                data.push((other.mat.row_indices()[k], i, other.mat.values()[k]));
            }
        }
        let triplet = Triplet {
            rows: self.mat.nrows().max(other.mat.nrows()),
            cols: self.mat.ncols().max(other.mat.ncols()),
            data,
        };
        SparseMatrix::from_triplet(triplet)
    }

    pub fn times_real(&self, s: f64) -> SparseMatrix {
        let mut data = Vec::new();
        for i in 0..self.mat.ncols() {
            let start = self.mat.col_ptrs()[i];
            let end = self.mat.col_ptrs()[i+1];
            for k in start..end {
                data.push((self.mat.row_indices()[k], i, self.mat.values()[k] * s));
            }
        }
        let triplet = Triplet {
            rows: self.mat.nrows(),
            cols: self.mat.ncols(),
            data,
        };
        SparseMatrix::from_triplet(triplet)
    }

    pub fn times_dense(&self, other: &DenseMatrix) -> DenseMatrix {
        let mut res = Mat::zeros(self.mat.nrows(), other.n_cols());
        for j in 0..other.n_cols() {
            for i in 0..self.mat.ncols() {
                let start = self.mat.col_ptrs()[i];
                let end = self.mat.col_ptrs()[i+1];
                let val_other = other.get(i, j);
                for k in start..end {
                    res[(self.mat.row_indices()[k], j)] += self.mat.values()[k] * val_other;
                }
            }
        }
        DenseMatrix { inner: res }
    }

    pub fn transpose(&self) -> SparseMatrix {
        let mut data = Vec::new();
        for i in 0..self.mat.ncols() {
            let start = self.mat.col_ptrs()[i];
            let end = self.mat.col_ptrs()[i+1];
            for k in start..end {
                data.push((i, self.mat.row_indices()[k], self.mat.values()[k]));
            }
        }
        let triplet = Triplet {
            rows: self.mat.ncols(),
            cols: self.mat.nrows(),
            data,
        };
        SparseMatrix::from_triplet(triplet)
    }

    pub fn invert_diagonal(&self) -> SparseMatrix {
        let mut data = Vec::new();
        for j in 0..self.mat.ncols() {
            let start = self.mat.col_ptrs()[j];
            let end = self.mat.col_ptrs()[j+1];
            for k in start..end {
                let i = self.mat.row_indices()[k];
                if i == j {
                    data.push((i, j, 1.0 / self.mat.values()[k]));
                }
            }
        }
        let triplet = Triplet {
            rows: self.mat.nrows(),
            cols: self.mat.ncols(),
            data,
        };
        SparseMatrix::from_triplet(triplet)
    }

    pub fn times_sparse(&self, other: &SparseMatrix) -> SparseMatrix {
        let mut t = Triplet::new(self.mat.nrows(), other.mat.ncols());
        for j in 0..other.mat.ncols() {
            let mut col_res = vec![0.0; self.mat.nrows()];
            let start_j = other.mat.col_ptrs()[j];
            let end_j = other.mat.col_ptrs()[j+1];
            for k_j in start_j..end_j {
                let mid_idx = other.mat.row_indices()[k_j];
                let val_other = other.mat.values()[k_j];
                
                let start_i = self.mat.col_ptrs()[mid_idx];
                let end_i = self.mat.col_ptrs()[mid_idx+1];
                for k_i in start_i..end_i {
                    col_res[self.mat.row_indices()[k_i]] += self.mat.values()[k_i] * val_other;
                }
            }
            for (i, &val) in col_res.iter().enumerate() {
                if val != 0.0 {
                    t.add_entry(val, i, j);
                }
            }
        }
        SparseMatrix::from_triplet(t)
    }

    pub fn identity(m: usize, n: usize) -> SparseMatrix {
        let mut t = Triplet::new(m, n);
        for i in 0..m.min(n) {
            t.add_entry(1.0, i, i);
        }
        SparseMatrix::from_triplet(t)
    }

    pub fn diag(v: &DenseMatrix) -> SparseMatrix {
        let n = v.n_rows();
        let mut t = Triplet::new(n, n);
        for i in 0..n {
            t.add_entry(v.get(i, 0), i, i);
        }
        SparseMatrix::from_triplet(t)
    }

    pub fn sub_matrix(&self, r0: usize, r1: usize, c0: usize, c1: usize) -> SparseMatrix {
        let mut t = Triplet::new(r1 - r0, c1 - c0);
        for j in c0..c1 {
            let start = self.mat.col_ptrs()[j];
            let end = self.mat.col_ptrs()[j+1];
            for k in start..end {
                let i = self.mat.row_indices()[k];
                if i >= r0 && i < r1 {
                    t.add_entry(self.mat.values()[k], i - r0, j - c0);
                }
            }
        }
        SparseMatrix::from_triplet(t)
    }

    pub fn increment_by(&mut self, other: &SparseMatrix) {
        *self = self.plus(other);
    }

    pub fn scale_by(&mut self, s: f64) {
        *self = self.times_real(s);
    }

    pub fn chol(&self) -> Cholesky {
        Cholesky { 
            mat: self.mat.to_dense()
        }
    }

    pub fn lu(&self) -> LU {
        LU {
            mat: self.mat.to_dense()
        }
    }

    pub fn qr(&self) -> QR {
        QR {
            mat: self.mat.to_dense()
        }
    }
}

pub struct ComplexSparseMatrix {
    pub mat: SparseColMat<usize, Complex64>,
}

impl ComplexSparseMatrix {
    pub fn from_triplet(triplet: ComplexTriplet) -> Self {
        let mat = SparseColMat::try_new_from_triplets(
            triplet.rows,
            triplet.cols,
            &triplet.data,
        ).unwrap();
        
        Self { mat }
    }

    pub fn n_rows(&self) -> usize { self.mat.nrows() }
    pub fn n_cols(&self) -> usize { self.mat.ncols() }

    pub fn scale_by(&mut self, s: Complex) {
        let mut data = Vec::new();
        let values = self.mat.values();
        for i in 0..self.mat.ncols() {
            let start = self.mat.col_ptrs()[i];
            let end = self.mat.col_ptrs()[i+1];
            for k in start..end {
                let val = Complex64::new(values.re[k], values.im[k]);
                data.push((self.mat.row_indices()[k], i, val * s.inner));
            }
        }
        let triplet = ComplexTriplet {
            rows: self.mat.nrows(),
            cols: self.mat.ncols(),
            data,
        };
        *self = ComplexSparseMatrix::from_triplet(triplet);
    }

    pub fn nnz(&self) -> usize {
        self.mat.compute_nnz()
    }

    pub fn frobenius_norm(&self) -> f64 {
        let mut sum_sq = 0.0;
        let values = self.mat.values();
        for i in 0..values.re.len() {
            sum_sq += values.re[i] * values.re[i] + values.im[i] * values.im[i];
        }
        sum_sq.sqrt()
    }

    pub fn plus(&self, other: &ComplexSparseMatrix) -> ComplexSparseMatrix {
        let mut data = Vec::new();
        let values = self.mat.values();
        for i in 0..self.mat.ncols() {
            let start = self.mat.col_ptrs()[i];
            let end = self.mat.col_ptrs()[i+1];
            for k in start..end {
                let val = Complex64::new(values.re[k], values.im[k]);
                data.push((self.mat.row_indices()[k], i, val));
            }
        }
        let values_other = other.mat.values();
        for i in 0..other.mat.ncols() {
            let start = other.mat.col_ptrs()[i];
            let end = other.mat.col_ptrs()[i+1];
            for k in start..end {
                let val = Complex64::new(values_other.re[k], values_other.im[k]);
                data.push((other.mat.row_indices()[k], i, val));
            }
        }
        let triplet = ComplexTriplet {
            rows: self.mat.nrows().max(other.mat.nrows()),
            cols: self.mat.ncols().max(other.mat.ncols()),
            data,
        };
        ComplexSparseMatrix::from_triplet(triplet)
    }

    pub fn minus(&self, other: &ComplexSparseMatrix) -> ComplexSparseMatrix {
        let mut data = Vec::new();
        let values = self.mat.values();
        for i in 0..self.mat.ncols() {
            let start = self.mat.col_ptrs()[i];
            let end = self.mat.col_ptrs()[i+1];
            for k in start..end {
                let val = Complex64::new(values.re[k], values.im[k]);
                data.push((self.mat.row_indices()[k], i, val));
            }
        }
        let values_other = other.mat.values();
        for i in 0..other.mat.ncols() {
            let start = other.mat.col_ptrs()[i];
            let end = other.mat.col_ptrs()[i+1];
            for k in start..end {
                let val = Complex64::new(values_other.re[k], values_other.im[k]);
                data.push((other.mat.row_indices()[k], i, -val));
            }
        }
        let triplet = ComplexTriplet {
            rows: self.mat.nrows().max(other.mat.nrows()),
            cols: self.mat.ncols().max(other.mat.ncols()),
            data,
        };
        ComplexSparseMatrix::from_triplet(triplet)
    }

    pub fn to_dense(&self) -> ComplexDenseMatrix {
        let mut res = ComplexDenseMatrix::zeros(self.n_rows(), self.n_cols());
        let values = self.mat.values();
        for j in 0..self.n_cols() {
            let start = self.mat.col_ptrs()[j];
            let end = self.mat.col_ptrs()[j+1];
            for k in start..end {
                let i = self.mat.row_indices()[k];
                res.set(Complex { inner: Complex64::new(values.re[k], values.im[k]) }, i, j);
            }
        }
        res
    }
}

pub struct Cholesky {
    mat: Mat<f64>,
}

impl Cholesky {
    pub fn solve_positive_definite(&self, rhs: &DenseMatrix) -> DenseMatrix {
        let llt = self.mat.as_ref().cholesky(faer::Side::Lower).unwrap();
        let mut res = rhs.inner.clone();
        llt.solve_in_place_with_conj_impl(res.as_mut(), faer::Conj::No);
        DenseMatrix { inner: res }
    }
}

pub struct LU {
    mat: Mat<f64>,
}

impl LU {
    pub fn solve_square(&self, rhs: &DenseMatrix) -> DenseMatrix {
        let lu = self.mat.as_ref().partial_piv_lu();
        let mut res = rhs.inner.clone();
        lu.solve_in_place_with_conj_impl(res.as_mut(), faer::Conj::No);
        DenseMatrix { inner: res }
    }
}

pub struct QR {
    mat: Mat<f64>,
}

impl QR {
    pub fn solve(&self, rhs: &DenseMatrix) -> DenseMatrix {
        let qr = self.mat.as_ref().col_piv_qr();
        let mut res = rhs.inner.clone();
        qr.solve_in_place_with_conj_impl(res.as_mut(), faer::Conj::No);
        DenseMatrix { inner: res }
    }
}
