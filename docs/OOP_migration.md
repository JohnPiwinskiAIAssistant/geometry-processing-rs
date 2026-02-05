# Object-Oriented Migration Guide

This document outlines a strategy for restructuring the Rust port to follow more object-oriented design principles while maintaining Rust's safety guarantees and performance characteristics.

## Current Architecture Analysis

### Current Design Pattern: Index-Based with Separate Storage

The current implementation uses a **data-oriented design** pattern:

```rust
// Current: Elements store indices, Mesh stores all data
pub struct Mesh {
    pub vertices: Vec<Vertex>,
    pub edges: Vec<Edge>,
    pub faces: Vec<Face>,
    pub halfedges: Vec<Halfedge>,
    // ...
}

pub struct Vertex {
    pub halfedge: Option<usize>,  // Index into mesh.halfedges
    pub index: usize,
}

pub struct Geometry<'a> {
    pub mesh: &'a Mesh,
    pub positions: Vec<Vector>,  // Separate from mesh
}
```

**Characteristics:**
- Elements are lightweight (just indices)
- All data centralized in `Mesh` and `Geometry`
- Methods require passing `&Mesh` or `&Geometry` references
- Cache-friendly for iteration
- Difficult to enforce invariants at compile time

### Issues with Current Design

1. **Tight Coupling**: Elements like `Vertex` require `Mesh` for most operations
2. **Scattered Behavior**: Geometric operations separated from mesh topology
3. **Weak Encapsulation**: Public fields expose internal structure
4. **Verbose API**: `vertex.adjacent_vertices(mesh, true)` instead of `vertex.adjacent_vertices()`

## Proposed OOP Architecture

### Design Goal: Smart Handles with Encapsulated Data

Migrate to a **handle-based OOP design** that provides:
- Encapsulated access to mesh data
- Methods that don't require passing mesh references
- Stronger type safety and invariant enforcement
- More intuitive API

### Core Concept: Smart Handles

```rust
// Proposed: Handles contain reference to mesh
pub struct VertexHandle<'a> {
    mesh: &'a Mesh,
    index: usize,
}

impl<'a> VertexHandle<'a> {
    // No need to pass mesh - it's already encapsulated!
    pub fn adjacent_vertices(&self) -> impl Iterator<Item = VertexHandle<'a>> + 'a {
        self.mesh.vertex_adjacent_vertices(self.index, true)
            .map(|idx| VertexHandle { mesh: self.mesh, index: idx })
    }
    
    pub fn position(&self, geometry: &Geometry) -> Vector {
        geometry.positions[self.index]
    }
    
    pub fn on_boundary(&self) -> bool {
        self.mesh.on_boundary(self.index)
    }
}
```

## Migration Strategy

### Phase 1: Introduce Handle Types

Create handle wrappers alongside existing types:

```rust
// In src/core/handles.rs
pub struct VertexHandle<'a> {
    mesh: &'a Mesh,
    index: usize,
}

pub struct EdgeHandle<'a> {
    mesh: &'a Mesh,
    index: usize,
}

pub struct FaceHandle<'a> {
    mesh: &'a Mesh,
    index: usize,
}

pub struct HalfedgeHandle<'a> {
    mesh: &'a Mesh,
    index: usize,
}
```

**Benefits:**
- Backward compatible (existing code still works)
- Gradual migration path
- Can be introduced incrementally

### Phase 2: Add Handle Methods to Mesh

```rust
impl Mesh {
    pub fn vertex(&self, index: usize) -> VertexHandle {
        VertexHandle { mesh: self, index }
    }
    
    pub fn vertices(&self) -> impl Iterator<Item = VertexHandle> + '_ {
        (0..self.vertices.len()).map(|i| self.vertex(i))
    }
    
    pub fn edge(&self, index: usize) -> EdgeHandle {
        EdgeHandle { mesh: self, index }
    }
    
    // ... similar for faces, halfedges
}
```

**Usage:**
```rust
// Old way
for v_idx in 0..mesh.vertices.len() {
    let neighbors = mesh.vertex_adjacent_vertices(v_idx, true);
}

// New way
for vertex in mesh.vertices() {
    let neighbors = vertex.adjacent_vertices();
}
```

### Phase 3: Integrate Geometry with Mesh

Currently, `Geometry` is separate from `Mesh`. Integrate them:

```rust
pub struct MeshWithGeometry {
    topology: Mesh,
    positions: Vec<Vector>,
}

impl MeshWithGeometry {
    pub fn vertex(&self, index: usize) -> GeometricVertexHandle {
        GeometricVertexHandle {
            mesh: &self.topology,
            geometry: self,
            index,
        }
    }
}

pub struct GeometricVertexHandle<'a> {
    mesh: &'a Mesh,
    geometry: &'a MeshWithGeometry,
    index: usize,
}

impl<'a> GeometricVertexHandle<'a> {
    pub fn position(&self) -> Vector {
        self.geometry.positions[self.index]
    }
    
    pub fn normal(&self) -> Vector {
        // Compute vertex normal using adjacent faces
        let mut n = Vector::new(0.0, 0.0, 0.0);
        for face in self.adjacent_faces() {
            n.increment_by(face.normal());
        }
        n.unit()
    }
    
    pub fn curvature(&self) -> f64 {
        // Compute using geometry methods
        self.geometry.scalar_gauss_curvature(self.index)
    }
}
```

### Phase 4: Encapsulate Fields

Make fields private and provide accessors:

```rust
pub struct VertexHandle<'a> {
    mesh: &'a Mesh,
    index: usize,  // Private
}

impl<'a> VertexHandle<'a> {
    pub fn index(&self) -> usize {
        self.index
    }
    
    pub fn degree(&self) -> usize {
        self.mesh.vertex_degree(self.index)
    }
    
    // Provide iterator-based access instead of exposing indices
    pub fn adjacent_vertices(&self) -> impl Iterator<Item = VertexHandle<'a>> + 'a {
        self.mesh.vertex_adjacent_vertices(self.index, true)
            .map(|idx| VertexHandle { mesh: self.mesh, index: idx })
    }
}
```

## Advanced OOP Patterns

### Pattern 1: Builder Pattern for Mesh Construction

```rust
pub struct MeshBuilder {
    vertices: Vec<Vector>,
    faces: Vec<[usize; 3]>,
}

impl MeshBuilder {
    pub fn new() -> Self {
        Self {
            vertices: Vec::new(),
            faces: Vec::new(),
        }
    }
    
    pub fn add_vertex(&mut self, pos: Vector) -> usize {
        let idx = self.vertices.len();
        self.vertices.push(pos);
        idx
    }
    
    pub fn add_face(&mut self, v0: usize, v1: usize, v2: usize) {
        self.faces.push([v0, v1, v2]);
    }
    
    pub fn build(self) -> Result<MeshWithGeometry, MeshError> {
        // Validate and build halfedge structure
        let mut mesh = Mesh::new();
        let soup = PolygonSoup {
            v: self.vertices,
            f: self.faces.into_iter().flatten().collect(),
        };
        
        if !mesh.build(&soup) {
            return Err(MeshError::NonManifold);
        }
        
        Ok(MeshWithGeometry {
            topology: mesh,
            positions: soup.v,
        })
    }
}
```

### Pattern 2: Visitor Pattern for Mesh Traversal

```rust
pub trait MeshVisitor {
    fn visit_vertex(&mut self, vertex: VertexHandle);
    fn visit_edge(&mut self, edge: EdgeHandle);
    fn visit_face(&mut self, face: FaceHandle);
}

impl Mesh {
    pub fn accept<V: MeshVisitor>(&self, visitor: &mut V) {
        for vertex in self.vertices() {
            visitor.visit_vertex(vertex);
        }
        for edge in self.edges() {
            visitor.visit_edge(edge);
        }
        for face in self.faces() {
            visitor.visit_face(face);
        }
    }
}

// Example: Compute mesh statistics
struct StatisticsVisitor {
    total_degree: usize,
    boundary_vertices: usize,
}

impl MeshVisitor for StatisticsVisitor {
    fn visit_vertex(&mut self, vertex: VertexHandle) {
        self.total_degree += vertex.degree();
        if vertex.on_boundary() {
            self.boundary_vertices += 1;
        }
    }
    
    fn visit_edge(&mut self, _edge: EdgeHandle) {}
    fn visit_face(&mut self, _face: FaceHandle) {}
}
```

### Pattern 3: Strategy Pattern for Geometric Operations

```rust
pub trait VertexNormalStrategy {
    fn compute(&self, vertex: &GeometricVertexHandle) -> Vector;
}

pub struct EquallyWeightedNormal;
impl VertexNormalStrategy for EquallyWeightedNormal {
    fn compute(&self, vertex: &GeometricVertexHandle) -> Vector {
        let mut n = Vector::new(0.0, 0.0, 0.0);
        for face in vertex.adjacent_faces() {
            n.increment_by(face.normal());
        }
        n.unit()
    }
}

pub struct AreaWeightedNormal;
impl VertexNormalStrategy for AreaWeightedNormal {
    fn compute(&self, vertex: &GeometricVertexHandle) -> Vector {
        let mut n = Vector::new(0.0, 0.0, 0.0);
        for face in vertex.adjacent_faces() {
            let normal = face.normal();
            let area = face.area();
            n.increment_by(normal.times(area));
        }
        n.unit()
    }
}

impl<'a> GeometricVertexHandle<'a> {
    pub fn normal_with_strategy<S: VertexNormalStrategy>(&self, strategy: &S) -> Vector {
        strategy.compute(self)
    }
}
```

### Pattern 4: Trait-Based Polymorphism

```rust
pub trait MeshElement {
    fn index(&self) -> usize;
    fn on_boundary(&self) -> bool;
}

pub trait GeometricElement: MeshElement {
    fn position(&self) -> Vector;
}

impl MeshElement for VertexHandle<'_> {
    fn index(&self) -> usize { self.index }
    fn on_boundary(&self) -> bool {
        self.mesh.on_boundary(self.index)
    }
}

impl GeometricElement for GeometricVertexHandle<'_> {
    fn position(&self) -> Vector {
        self.geometry.positions[self.index]
    }
}

// Generic function that works with any mesh element
fn print_element_info<E: MeshElement>(element: &E) {
    println!("Element {} is {}on boundary",
        element.index(),
        if element.on_boundary() { "" } else { "not " }
    );
}
```

## Migration Checklist

### Step 1: Create Handle Types
- [ ] Create `src/core/handles.rs`
- [ ] Define `VertexHandle`, `EdgeHandle`, `FaceHandle`, `HalfedgeHandle`
- [ ] Implement basic accessors

### Step 2: Add Handle Constructors to Mesh
- [ ] Add `mesh.vertex(index)` method
- [ ] Add `mesh.vertices()` iterator
- [ ] Similar for edges, faces, halfedges

### Step 3: Migrate Adjacency Methods
- [ ] Move `vertex.adjacent_vertices()` to return handles
- [ ] Move `edge.adjacent_faces()` to return handles
- [ ] Update all adjacency methods

### Step 4: Integrate Geometry
- [ ] Create `MeshWithGeometry` struct
- [ ] Create `GeometricVertexHandle`, etc.
- [ ] Move geometric methods to handles

### Step 5: Update Projects
- [ ] Refactor `Parameterization` to use handles
- [ ] Refactor `PoissonProblem` to use handles
- [ ] Update all other projects

### Step 6: Encapsulation
- [ ] Make handle fields private
- [ ] Provide controlled access through methods
- [ ] Add validation where appropriate

## Trade-offs and Considerations

### Advantages of OOP Migration

1. **Better API**: More intuitive, less verbose
2. **Encapsulation**: Hide implementation details
3. **Type Safety**: Handles prevent invalid indices
4. **Discoverability**: IDE autocomplete works better
5. **Maintainability**: Changes localized to handle implementations

### Disadvantages

1. **Performance**: Extra indirection (minimal in practice)
2. **Lifetime Complexity**: Handles tied to mesh lifetime
3. **Migration Effort**: Significant refactoring required
4. **Memory**: Handles larger than raw indices (but stack-allocated)

### Rust-Specific Considerations

**Borrowing Rules:**
```rust
// This won't compile - can't mutate mesh while handles exist
let vertex = mesh.vertex(0);
mesh.add_vertex(Vector::new(1.0, 0.0, 0.0));  // Error!
```

**Solution:** Use interior mutability or separate mutable operations:
```rust
impl Mesh {
    pub fn modify<F>(&mut self, f: F)
    where F: FnOnce(&mut Mesh)
    {
        f(self);
    }
}
```

## Conclusion

This migration to OOP patterns will make the codebase more maintainable and user-friendly while preserving Rust's safety guarantees. The handle-based approach is idiomatic Rust and aligns well with how modern Rust libraries (like `petgraph`) structure their APIs.

The migration can be done incrementally, allowing the old and new APIs to coexist during the transition period.
