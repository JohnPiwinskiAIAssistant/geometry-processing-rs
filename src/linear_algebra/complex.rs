use num_complex::Complex64;
use std::ops::{Add, Sub, Mul, Div, Neg};

#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Complex {
    pub inner: Complex64,
}

impl Complex {
    pub fn new(re: f64, im: f64) -> Self {
        Self {
            inner: Complex64::new(re, im),
        }
    }

    pub fn re(&self) -> f64 {
        self.inner.re
    }

    pub fn im(&self) -> f64 {
        self.inner.im
    }

    pub fn arg(&self) -> f64 {
        self.inner.arg()
    }

    pub fn norm(&self) -> f64 {
        self.inner.norm()
    }

    pub fn norm2(&self) -> f64 {
        self.inner.norm_sqr()
    }

    pub fn conjugate(&self) -> Self {
        Self {
            inner: self.inner.conj(),
        }
    }

    pub fn inverse(&self) -> Self {
        Self {
            inner: self.inner.inv(),
        }
    }

    pub fn polar(&self) -> Self {
        let (r, theta) = self.inner.to_polar();
        Self::new(theta.cos() * r, theta.sin() * r)
    }

    pub fn exp(&self) -> Self {
        Self {
            inner: self.inner.exp(),
        }
    }

    pub fn plus(&self, v: Complex) -> Self {
        *self + v
    }

    pub fn minus(&self, v: Complex) -> Self {
        *self - v
    }

    pub fn times_real(&self, s: f64) -> Self {
        Self {
            inner: self.inner * s,
        }
    }

    pub fn over_real(&self, s: f64) -> Self {
        Self {
            inner: self.inner / s,
        }
    }

    pub fn times_complex(&self, v: Complex) -> Self {
        Self {
            inner: self.inner * v.inner,
        }
    }

    pub fn over_complex(&self, v: Complex) -> Self {
        Self {
            inner: self.inner / v.inner,
        }
    }
}

// Trait implementations
impl Add for Complex {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        Self { inner: self.inner + rhs.inner }
    }
}

impl Sub for Complex {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self::Output {
        Self { inner: self.inner - rhs.inner }
    }
}

impl Mul for Complex {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        Self { inner: self.inner * rhs.inner }
    }
}

impl Div for Complex {
    type Output = Self;
    fn div(self, rhs: Self) -> Self::Output {
        Self { inner: self.inner / rhs.inner }
    }
}

impl Neg for Complex {
    type Output = Self;
    fn neg(self) -> Self::Output {
        Self { inner: -self.inner }
    }
}
