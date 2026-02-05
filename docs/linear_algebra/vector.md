# Vector

The `Vector` struct represents an element of Euclidean 3-space $(x, y, z)$, along with all the usual vector space operations.

```rust
pub struct Vector {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}
```

This is a fundamental type used throughout the geometry processing framework for representing positions, normals, and other geometric quantities.

## Constructor

```rust
pub fn new(x: f64, y: f64, z: f64) -> Self
```

Creates a new vector with the specified components. Corresponds to the JavaScript constructor `constructor(x = 0, y = 0, z = 0)`.

## Norms and Normalization

```rust
pub fn norm(&self) -> f64
```

Computes the Euclidean length $\sqrt{x^2 + y^2 + z^2}$ of this vector.

---

```rust
pub fn norm2(&self) -> f64
```

Computes the squared Euclidean length $x^2 + y^2 + z^2$. This is more efficient when only relative magnitudes are needed.

---

```rust
pub fn normalize(&mut self)
```

Divides this vector by its Euclidean length in-place, making it a unit vector.

---

```rust
pub fn unit(&self) -> Self
```

Returns a normalized copy of this vector without modifying the original.

---

```rust
pub fn is_valid(&self) -> bool
```

Checks whether all components are finite (not NaN or infinite). Useful for debugging geometric algorithms.

## In-Place Operations

These methods modify the vector in place:

```rust
pub fn increment_by(&mut self, v: Vector)  // u += v
pub fn decrement_by(&mut self, v: Vector)  // u -= v
pub fn scale_by(&mut self, s: f64)         // u *= s
pub fn divide_by(&mut self, s: f64)        // u /= s
```

## Arithmetic Operations

These methods return new vectors:

```rust
pub fn plus(&self, v: Vector) -> Self      // u + v
pub fn minus(&self, v: Vector) -> Self     // u - v
pub fn times(&self, s: f64) -> Self        // u * s
pub fn over(&self, s: f64) -> Self         // u / s
pub fn negated(&self) -> Self              // -u
```

## Vector Products

```rust
pub fn dot(&self, v: Vector) -> f64
```

Computes the dot product $u \cdot v = u_x v_x + u_y v_y + u_z v_z$.

---

```rust
pub fn cross(&self, v: Vector) -> Self
```

Computes the cross product $u \times v$, which is perpendicular to both $u$ and $v$ with magnitude equal to the area of the parallelogram they span.

## Trait Implementations

For idiomatic Rust usage, `Vector` implements:
- `Add`, `Sub`: `+`, `-` operators
- `Mul<f64>`, `Div<f64>`: `*`, `/` with scalars
- `AddAssign`, `SubAssign`, `MulAssign`, `DivAssign`: `+=`, `-=`, `*=`, `/=`
- `Neg`: Unary `-`
- `Clone`, `Copy`, `Debug`, `PartialEq`, `Default`
