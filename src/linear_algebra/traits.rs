use faer::Mat;
use num_complex::Complex64;
use std::ops::{Add, Sub, Mul, Div, Neg};

/// Trait for scalars (f64 or Complex64)
pub trait Scalar: 
    Clone + Copy + 
    Add<Output = Self> + Sub<Output = Self> + Mul<Output = Self> + Div<Output = Self> + 
    Neg<Output = Self> +
    PartialEq +
    Default +
    faer::Entity +
    'static
{
    fn zero() -> Self;
    fn one() -> Self;
    fn from_f64(v: f64) -> Self;
    fn abs(&self) -> f64;
    fn sqrt(&self) -> Self;
}

impl Scalar for f64 {
    fn zero() -> Self { 0.0 }
    fn one() -> Self { 1.0 }
    fn from_f64(v: f64) -> Self { v }
    fn abs(&self) -> f64 { f64::abs(*self) }
    fn sqrt(&self) -> Self { f64::sqrt(*self) }
}

impl Scalar for Complex64 {
    fn zero() -> Self { Complex64::new(0.0, 0.0) }
    fn one() -> Self { Complex64::new(1.0, 0.0) }
    fn from_f64(v: f64) -> Self { Complex64::new(v, 0.0) }
    fn abs(&self) -> f64 { Complex64::norm(*self) }
    fn sqrt(&self) -> Self { Complex64::sqrt(*self) }
}

/// Generic matrix operations for dense matrices/vectors
pub trait DenseMatrixOps<S: Scalar>: Sized {
    fn nrows(&self) -> usize;
    fn ncols(&self) -> usize;
    fn get(&self, row: usize, col: usize) -> S;
    fn set(&mut self, row: usize, col: usize, val: S);
    
    fn norm_sq(&self) -> f64;
    fn norm(&self) -> f64 { self.norm_sq().sqrt() }
    
    fn dot(&self, other: &Self) -> S;
    fn transpose(&self) -> Self;
}

impl<S: Scalar> DenseMatrixOps<S> for Mat<S> {
    fn nrows(&self) -> usize { self.nrows() }
    fn ncols(&self) -> usize { self.ncols() }
    
    fn get(&self, row: usize, col: usize) -> S { self.read(row, col) }
    fn set(&mut self, row: usize, col: usize, val: S) { self.write(row, col, val); }
    
    fn norm_sq(&self) -> f64 {
        let mut sum = 0.0;
        for i in 0..self.nrows() {
            for j in 0..self.ncols() {
                let val = self.read(i, j);
                sum += val.abs() * val.abs();
            }
        }
        sum
    }
    
    fn dot(&self, other: &Self) -> S {
        assert_eq!(self.nrows(), other.nrows());
        assert_eq!(self.ncols(), other.ncols());
        let mut sum = S::zero();
        for i in 0..self.nrows() {
            for j in 0..self.ncols() {
                sum = sum + self.read(i, j) * other.read(i, j);
            }
        }
        sum
    }
    
    fn transpose(&self) -> Self {
        Mat::from_fn(self.ncols(), self.nrows(), |i, j| self.read(j, i))
    }
}

/// Operations specific to 3D geometric vectors
pub trait Vector3Ops<S: Scalar>: DenseMatrixOps<S> {
    fn cross(&self, other: &Self) -> Self;
    fn unit(&self) -> Self;
}

impl<S: Scalar> Vector3Ops<S> for Mat<S> {
    fn cross(&self, other: &Self) -> Self {
        assert_eq!(self.nrows(), 3);
        assert_eq!(self.ncols(), 1);
        assert_eq!(other.nrows(), 3);
        assert_eq!(other.ncols(), 1);
        
        let a1 = self.read(0, 0);
        let a2 = self.read(1, 0);
        let a3 = self.read(2, 0);
        
        let b1 = other.read(0, 0);
        let b2 = other.read(1, 0);
        let b3 = other.read(2, 0);
        
        faer::mat![
            [a2 * b3 - a3 * b2],
            [a3 * b1 - a1 * b3],
            [a1 * b2 - a2 * b1]
        ]
    }
    
    fn unit(&self) -> Self {
        let n = self.norm();
        if n == 0.0 {
            self.clone()
        } else {
            let inv_n = S::from_f64(1.0 / n);
            let mut res = self.clone();
            for i in 0..res.nrows() {
                for j in 0..res.ncols() {
                    let val = res.read(i, j);
                    res.write(i, j, val * inv_n);
                }
            }
            res
        }
    }
}

/// Operations for sparse matrices
pub trait SparseOps<S: Scalar>: Sized {
    fn from_triplets(rows: usize, cols: usize, triplets: &[(usize, usize, S)]) -> Self;
    fn invert_diagonal(&self) -> Self;
    fn sub_matrix(&self, r0: usize, r1: usize, c0: usize, c1: usize) -> Self;
    fn get_val(&self, row: usize, col: usize) -> S;
    fn frobenius_norm(&self) -> f64;
    fn scale(&self, s: S) -> Self;
    fn compute_nnz(&self) -> usize;
}

/// Generic linear solver trait
pub trait LinearSolver<S: Scalar> {
    fn solve(&self, rhs: &Mat<S>) -> Mat<S>;
}
