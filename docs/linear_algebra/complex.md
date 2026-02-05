# Complex.rs Mapping

| Rust Function | JS Function | JS Line Number | Description | Notes |
|---|---|---|---|---|
| `Complex::new` | `Complex` (constructor) | 12 | Initializes a complex number (real, imag). | Internal storage uses `num_complex::Complex64` |
| `Complex::re` | `Complex.re` | 21 | Returns the real part. | |
| `Complex::im` | `Complex.im` | 30 | Returns the imaginary part. | |
| `Complex::arg` | `Complex.arg` | 39 | Computes the argument (phase). | |
| `Complex::norm` | `Complex.norm` | 48 | Computes the modulus (magnitude). | |
| `Complex::norm2` | `Complex.norm2` | 57 | Computes the squared modulus. | |
| `Complex::conjugate` | `Complex.conjugate` | 66 | Returns the complex conjugate. | |
| `Complex::inverse` | `Complex.inverse` | 75 | Returns the multiplicative inverse. | |
| `Complex::polar` | (No direct match) | | Constructs from polar coords (r, theta). | Rust implementation detail (or helper) |
| `Complex::exp` | `Complex.exp` | 84 | Computes the complex exponential. | |
| `Complex::plus` | `Complex.plus` | 94 | Adds another complex number. | |
| `Complex::minus` | `Complex.minus` | 105 | Subtracts another complex number. | |
| `Complex::times_real` | `Complex.times` (scalar overload) | 116 | Multiplies by a real scalar. | JS `times` handles both scalar and complex |
| `Complex::over_real` | `Complex.over` (scalar overload) | 127 | Divides by a real scalar. | JS `over` handles both scalar and complex |
| `Complex::times_complex` | `Complex.times` (complex overload) | 116 | Multiplies by a complex number. | |
| `Complex::over_complex` | `Complex.over` (complex overload) | 127 | Divides by a complex number. | |
| (Trait Impl) `Add` | | | Operator overload for `+`. | |
| (Trait Impl) `Sub` | | | Operator overload for `-`. | |
| (Trait Impl) `Mul` | | | Operator overload for `*`. | |
| (Trait Impl) `Div` | | | Operator overload for `/`. | |
| (Trait Impl) `Neg` | | | Operator overload for `-` (unary). | |
