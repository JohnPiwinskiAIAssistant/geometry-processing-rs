# Corner.rs Mapping

| Rust Function | JS Function | JS Line Number | Description | Notes |
|---|---|---|---|---|
| `Corner::new` | `Corner.constructor` | 10 | Initializes a new Corner instance. | |
| `Corner::vertex` | `Corner.vertex` | 20 | Gets the vertex this corner lies on. | Takes `&Mesh`; JS is a getter |
| `Corner::face` | `Corner.face` | 29 | Gets the face this corner is contained in. | Takes `&Mesh`; JS is a getter |
| `Corner::next` | `Corner.next` | 38 | Gets the next corner (CCW) in this face. | Takes `&Mesh`; JS is a getter |
| `Corner::prev` | `Corner.prev` | 47 | Gets the previous corner (CCW) in this face. | Takes `&Mesh`; JS is a getter |
| (Not implemented) | `Corner.toString` | 57 | Returns a string representation (index). | Implemented via `Debug` trait |
