use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Vector {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Vector {
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    pub fn norm(&self) -> f64 {
        self.norm2().sqrt()
    }

    pub fn norm2(&self) -> f64 {
        self.dot(*self)
    }

    pub fn normalize(&mut self) {
        let n = self.norm();
        self.x /= n;
        self.y /= n;
        self.z /= n;
    }

    pub fn unit(&self) -> Self {
        let n = self.norm();
        Self {
            x: self.x / n,
            y: self.y / n,
            z: self.z / n,
        }
    }

    pub fn is_valid(&self) -> bool {
        self.x.is_finite() && self.y.is_finite() && self.z.is_finite()
    }

    pub fn increment_by(&mut self, v: Vector) {
        *self += v;
    }

    pub fn decrement_by(&mut self, v: Vector) {
        *self -= v;
    }

    pub fn scale_by(&mut self, s: f64) {
        *self *= s;
    }

    pub fn divide_by(&mut self, s: f64) {
        *self /= s;
    }

    pub fn plus(&self, v: Vector) -> Self {
        *self + v
    }

    pub fn minus(&self, v: Vector) -> Self {
        *self - v
    }

    pub fn times(&self, s: f64) -> Self {
        *self * s
    }

    pub fn over(&self, s: f64) -> Self {
        *self / s
    }

    pub fn negated(&self) -> Self {
        -*self
    }

    pub fn dot(&self, v: Vector) -> f64 {
        self.x * v.x + self.y * v.y + self.z * v.z
    }

    pub fn cross(&self, v: Vector) -> Self {
        Self {
            x: self.y * v.z - self.z * v.y,
            y: self.z * v.x - self.x * v.z,
            z: self.x * v.y - self.y * v.x,
        }
    }
}

// Trait implementations for idiomatic Rust usage
impl Add for Vector {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        Self::new(self.x + rhs.x, self.y + rhs.y, self.z + rhs.z)
    }
}

impl AddAssign for Vector {
    fn add_assign(&mut self, rhs: Self) {
        self.x += rhs.x;
        self.y += rhs.y;
        self.z += rhs.z;
    }
}

impl Sub for Vector {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self::Output {
        Self::new(self.x - rhs.x, self.y - rhs.y, self.z - rhs.z)
    }
}

impl SubAssign for Vector {
    fn sub_assign(&mut self, rhs: Self) {
        self.x -= rhs.x;
        self.y -= rhs.y;
        self.z -= rhs.z;
    }
}

impl Mul<f64> for Vector {
    type Output = Self;
    fn mul(self, rhs: f64) -> Self::Output {
        Self::new(self.x * rhs, self.y * rhs, self.z * rhs)
    }
}

impl MulAssign<f64> for Vector {
    fn mul_assign(&mut self, rhs: f64) {
        self.x *= rhs;
        self.y *= rhs;
        self.z *= rhs;
    }
}

impl Div<f64> for Vector {
    type Output = Self;
    fn div(self, rhs: f64) -> Self::Output {
        Self::new(self.x / rhs, self.y / rhs, self.z / rhs)
    }
}

impl DivAssign<f64> for Vector {
    fn div_assign(&mut self, rhs: f64) {
        self.x /= rhs;
        self.y /= rhs;
        self.z /= rhs;
    }
}

impl Neg for Vector {
    type Output = Self;
    fn neg(self) -> Self::Output {
        Self::new(-self.x, -self.y, -self.z)
    }
}
