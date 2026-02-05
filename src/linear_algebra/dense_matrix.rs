use faer::Mat;
use num_complex::Complex64;
use crate::linear_algebra::Complex;

#[derive(Debug, Clone, PartialEq)]
pub struct DenseMatrix {
    pub inner: Mat<f64>,
}

impl DenseMatrix {
    pub fn from_inner(inner: Mat<f64>) -> Self {
        Self { inner }
    }

    pub fn zeros(m: usize, n: usize) -> Self {
        Self {
            inner: Mat::zeros(m, n),
        }
    }

    pub fn identity(m: usize, n: usize) -> Self {
        let mut mat = Mat::zeros(m, n);
        for i in 0..m.min(n) {
            mat[(i, i)] = 1.0;
        }
        Self { inner: mat }
    }

    pub fn ones(m: usize, n: usize) -> Self {
        Self {
            inner: Mat::from_fn(m, n, |_, _| 1.0),
        }
    }

    pub fn constant(x: f64, m: usize, n: usize) -> Self {
        Self {
            inner: Mat::from_fn(m, n, |_, _| x),
        }
    }

    pub fn random(m: usize, n: usize) -> Self {
        Self {
            inner: Mat::from_fn(m, n, |_, _| rand::random::<f64>()),
        }
    }

    pub fn transpose(&self) -> Self {
        Self {
            inner: self.inner.transpose().to_owned(),
        }
    }

    pub fn n_rows(&self) -> usize {
        self.inner.nrows()
    }

    pub fn n_cols(&self) -> usize {
        self.inner.ncols()
    }

    pub fn norm(&self, n: usize) -> f64 {
        match n {
            1 => {
                let mut max_col_sum = 0.0;
                for j in 0..self.n_cols() {
                    let mut col_sum = 0.0;
                    for i in 0..self.n_rows() {
                        col_sum += self.inner[(i, j)].abs();
                    }
                    if col_sum > max_col_sum { max_col_sum = col_sum; }
                }
                max_col_sum
            }
            0 => {
                let mut max_row_sum = 0.0;
                for i in 0..self.n_rows() {
                    let mut row_sum = 0.0;
                    for j in 0..self.n_cols() {
                        row_sum += self.inner[(i, j)].abs();
                    }
                    if row_sum > max_row_sum { max_row_sum = row_sum; }
                }
                max_row_sum
            }
            _ => {
                let mut sum_sq = 0.0;
                for j in 0..self.n_cols() {
                    for i in 0..self.n_rows() {
                        sum_sq += self.inner[(i, j)] * self.inner[(i, j)];
                    }
                }
                sum_sq.sqrt()
            }
        }
    }

    pub fn rank(&self) -> usize {
        self.inner.ncols() // stub
    }

    pub fn sum(&self) -> f64 {
        let mut s = 0.0;
        for j in 0..self.n_cols() {
            for i in 0..self.n_rows() {
                s += self.inner[(i, j)];
            }
        }
        s
    }

    pub fn sub_matrix(&self, r0: usize, r1: usize, c0: usize, c1: usize) -> Self {
        Self {
            inner: self.inner.as_ref().submatrix(r0, c0, r1 - r0, c1 - c0).to_owned(),
        }
    }

    pub fn increment_by(&mut self, other: &DenseMatrix) {
        self.inner += &other.inner;
    }

    pub fn decrement_by(&mut self, other: &DenseMatrix) {
        self.inner -= &other.inner;
    }

    pub fn scale_by(&mut self, s: f64) {
        self.inner *= s;
    }

    pub fn plus(&self, other: &DenseMatrix) -> Self {
        Self {
            inner: &self.inner + &other.inner,
        }
    }

    pub fn minus(&self, other: &DenseMatrix) -> Self {
        Self {
            inner: &self.inner - &other.inner,
        }
    }

    pub fn times_real(&self, s: f64) -> Self {
        Self {
            inner: &self.inner * s,
        }
    }

    pub fn times_dense(&self, other: &DenseMatrix) -> Self {
        Self {
            inner: &self.inner * &other.inner,
        }
    }

    pub fn negated(&self) -> Self {
        Self {
            inner: -&self.inner,
        }
    }

    pub fn get(&self, i: usize, j: usize) -> f64 {
        self.inner[(i, j)]
    }

    pub fn set(&mut self, x: f64, i: usize, j: usize) {
        self.inner[(i, j)] = x;
    }

    pub fn hcat(&self, other: &DenseMatrix) -> Self {
        let m = self.n_rows();
        let n1 = self.n_cols();
        let n2 = other.n_cols();
        let mut res = Mat::zeros(m, n1 + n2);
        for j in 0..n1 {
            for i in 0..m {
                res[(i, j)] = self.inner[(i, j)];
            }
        }
        for j in 0..n2 {
            for i in 0..m {
                res[(i, n1 + j)] = other.inner[(i, j)];
            }
        }
        Self { inner: res }
    }

    pub fn vcat(&self, other: &DenseMatrix) -> Self {
        let m1 = self.n_rows();
        let m2 = other.n_rows();
        let n = self.n_cols();
        let mut res = Mat::zeros(m1 + m2, n);
        for j in 0..n {
            for i in 0..m1 {
                res[(i, j)] = self.inner[(i, j)];
            }
            for i in 0..m2 {
                res[(m1 + i, j)] = other.inner[(i, j)];
            }
        }
        Self { inner: res }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct ComplexDenseMatrix {
    pub rows: usize,
    pub cols: usize,
    pub data: Vec<Complex64>,
}

impl ComplexDenseMatrix {
    pub fn zeros(m: usize, n: usize) -> Self {
        Self {
            rows: m,
            cols: n,
            data: vec![Complex64::new(0.0, 0.0); m * n],
        }
    }

    pub fn ones(m: usize, n: usize) -> Self {
        Self {
            rows: m,
            cols: n,
            data: vec![Complex64::new(1.0, 0.0); m * n],
        }
    }

    pub fn random(m: usize, n: usize) -> Self {
        let mut data = Vec::with_capacity(m * n);
        for _ in 0..m * n {
            data.push(Complex64::new(rand::random::<f64>(), rand::random::<f64>()));
        }
        Self {
            rows: m,
            cols: n,
            data,
        }
    }

    pub fn n_rows(&self) -> usize { self.rows }
    pub fn n_cols(&self) -> usize { self.cols }

    pub fn transpose(&self) -> Self {
        let mut res = Self::zeros(self.cols, self.rows);
        for j in 0..self.cols {
            for i in 0..self.rows {
                res.set(self.get(i, j), j, i);
            }
        }
        res
    }

    pub fn conjugate(&self) -> Self {
        let mut res = self.clone();
        for x in &mut res.data {
            *x = x.conj();
        }
        res
    }

    pub fn get(&self, i: usize, j: usize) -> Complex {
        Complex { inner: self.data[j * self.rows + i] }
    }

    pub fn set(&mut self, x: Complex, i: usize, j: usize) {
        self.data[j * self.rows + i] = x.inner;
    }

    pub fn sum(&self) -> Complex {
        let s = self.data.iter().sum();
        Complex { inner: s }
    }

    pub fn norm(&self, _n: usize) -> f64 {
        let sum_sq: f64 = self.data.iter().map(|x| x.norm_sqr()).sum();
        sum_sq.sqrt()
    }

    pub fn decrement_by(&mut self, other: &ComplexDenseMatrix) {
        for (i, x) in self.data.iter_mut().enumerate() {
            *x -= other.data[i];
        }
    }

    pub fn scale_by(&mut self, s: Complex) {
        for x in &mut self.data {
            *x *= s.inner;
        }
    }

    pub fn times_dense(&self, other: &ComplexDenseMatrix) -> Self {
        // Simple O(N^3) multiplication
        let mut res = Self::zeros(self.rows, other.cols);
        for j in 0..other.cols {
            for k in 0..self.cols {
                let val_other = other.get(k, j);
                for i in 0..self.rows {
                    let mut val = res.get(i, j);
                    val.inner += self.get(i, k).inner * val_other.inner;
                    res.set(val, i, j);
                }
            }
        }
        res
    }

    pub fn minus(&self, other: &ComplexDenseMatrix) -> Self {
        let mut res = self.clone();
        res.decrement_by(other);
        res
    }

    pub fn times_complex(&self, s: Complex) -> Self {
        let mut res = self.clone();
        res.scale_by(s);
        res
    }
}
