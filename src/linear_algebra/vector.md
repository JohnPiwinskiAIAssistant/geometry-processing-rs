# Vector.rs Mapping

| Rust Function | JS Function | JS Line Number | Description | Notes |
|---|---|---|---|---|
| `Vector::new` | `Vector.constructor` | 10 | Initializes a new Vector with x, y, z components. | |
| `Vector::norm` | `Vector.norm` | 22 | Computes the Euclidean norm (length). | |
| `Vector::norm2` | `Vector.norm2` | 31 | Computes the squared Euclidean norm. | |
| `Vector::normalize` | `Vector.normalize` | 39 | Normalizes the vector in place. | |
| `Vector::unit` | `Vector.unit` | 51 | Returns a new normalized unit vector. | |
| `Vector::is_valid` | `Vector.isValid` | 61 | Checks if all components are finite numbers. | |
| `Vector::increment_by` | `Vector.incrementBy` | 71 | Adds another vector to this one in place. | |
| `Vector::decrement_by` | `Vector.decrementBy` | 81 | Subtracts another vector from this one in place. | |
| `Vector::scale_by` | `Vector.scaleBy` | 91 | Scales this vector by a scalar in place. | |
| `Vector::divide_by` | `Vector.divideBy` | 101 | Divides this vector by a scalar in place. | |
| `Vector::plus` | `Vector.plus` | 114 | Returns sum of this vector and another. | |
| `Vector::minus` | `Vector.minus` | 125 | Returns difference of this vector and another. | |
| `Vector::times` | `Vector.times` | 136 | Returns this vector scaled by a scalar. | |
| `Vector::over` | `Vector.over` | 147 | Returns this vector divided by a scalar. | |
| `Vector::negated` | `Vector.negated` | 158 | Returns the negation of this vector. | |
| `Vector::dot` | `Vector.dot` | 169 | Computes dot product with another vector. | |
| `Vector::cross` | `Vector.cross` | 180 | Computes cross product with another vector. | |
| (Trait Impl) `Add` | | | Operator overload for `+`. | Rust idiomatic |
| (Trait Impl) `Sub` | | | Operator overload for `-`. | Rust idiomatic |
| (Trait Impl) `Mul` | | | Operator overload for `*`. | Rust idiomatic |
| (Trait Impl) `Div` | | | Operator overload for `/`. | Rust idiomatic |
| (Trait Impl) `Neg` | | | Operator overload for `-` (unary). | Rust idiomatic |
| (Trait Impl) `AddAssign` | | | Operator overload for `+=`. | Rust idiomatic |
| (Trait Impl) `SubAssign` | | | Operator overload for `-=`. | Rust idiomatic |
| (Trait Impl) `MulAssign` | | | Operator overload for `*=`. | Rust idiomatic |
| (Trait Impl) `DivAssign` | | | Operator overload for `/=`. | Rust idiomatic |
