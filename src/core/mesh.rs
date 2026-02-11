pub use crate::core::vertex::Vertex;
pub use crate::core::edge::Edge;
pub use crate::core::face::Face;
pub use crate::core::halfedge::Halfedge;
pub use crate::core::corner::Corner;
pub use crate::linear_algebra::{Vector, Complex, SparseMatrix, Triplet};
use std::collections::HashMap;

pub trait MeshBackend {
    fn num_vertices(&self) -> usize;
    fn num_edges(&self) -> usize;
    fn num_faces(&self) -> usize;
    fn num_halfedges(&self) -> usize;
    fn num_corners(&self) -> usize;
    fn num_boundaries(&self) -> usize;
    fn num_generators(&self) -> usize;

    fn vertex_halfedge(&self, v_idx: usize) -> Option<usize>;
    fn edge_halfedge(&self, e_idx: usize) -> Option<usize>;
    fn face_halfedge(&self, f_idx: usize) -> Option<usize>;
    fn corner_halfedge(&self, c_idx: usize) -> Option<usize>;
    fn boundary_halfedge(&self, b_idx: usize) -> Option<usize>;
    fn generator(&self, g_idx: usize) -> &[usize];

    fn halfedge_vertex(&self, h_idx: usize) -> Option<usize>;
    fn halfedge_edge(&self, h_idx: usize) -> Option<usize>;
    fn halfedge_face(&self, h_idx: usize) -> Option<usize>;
    fn halfedge_corner(&self, h_idx: usize) -> Option<usize>;
    fn halfedge_next(&self, h_idx: usize) -> Option<usize>;
    fn halfedge_prev(&self, h_idx: usize) -> Option<usize>;
    fn halfedge_twin(&self, h_idx: usize) -> Option<usize>;
    fn halfedge_on_boundary(&self, h_idx: usize) -> bool;

    fn add_generator(&mut self, generator: Vec<usize>);
}

impl<B: MeshBackend> MeshBackend for &B {
    fn num_vertices(&self) -> usize { (**self).num_vertices() }
    fn num_edges(&self) -> usize { (**self).num_edges() }
    fn num_faces(&self) -> usize { (**self).num_faces() }
    fn num_halfedges(&self) -> usize { (**self).num_halfedges() }
    fn num_corners(&self) -> usize { (**self).num_corners() }
    fn num_boundaries(&self) -> usize { (**self).num_boundaries() }
    fn num_generators(&self) -> usize { (**self).num_generators() }

    fn vertex_halfedge(&self, v_idx: usize) -> Option<usize> { (**self).vertex_halfedge(v_idx) }
    fn edge_halfedge(&self, e_idx: usize) -> Option<usize> { (**self).edge_halfedge(e_idx) }
    fn face_halfedge(&self, f_idx: usize) -> Option<usize> { (**self).face_halfedge(f_idx) }
    fn corner_halfedge(&self, c_idx: usize) -> Option<usize> { (**self).corner_halfedge(c_idx) }
    fn boundary_halfedge(&self, b_idx: usize) -> Option<usize> { (**self).boundary_halfedge(b_idx) }
    fn generator(&self, g_idx: usize) -> &[usize] { (**self).generator(g_idx) }

    fn halfedge_vertex(&self, h_idx: usize) -> Option<usize> { (**self).halfedge_vertex(h_idx) }
    fn halfedge_edge(&self, h_idx: usize) -> Option<usize> { (**self).halfedge_edge(h_idx) }
    fn halfedge_face(&self, h_idx: usize) -> Option<usize> { (**self).halfedge_face(h_idx) }
    fn halfedge_corner(&self, h_idx: usize) -> Option<usize> { (**self).halfedge_corner(h_idx) }
    fn halfedge_next(&self, h_idx: usize) -> Option<usize> { (**self).halfedge_next(h_idx) }
    fn halfedge_prev(&self, h_idx: usize) -> Option<usize> { (**self).halfedge_prev(h_idx) }
    fn halfedge_twin(&self, h_idx: usize) -> Option<usize> { (**self).halfedge_twin(h_idx) }
    fn halfedge_on_boundary(&self, h_idx: usize) -> bool { (**self).halfedge_on_boundary(h_idx) }

    fn add_generator(&mut self, _generator: Vec<usize>) {
        panic!("Cannot add generator to a MeshBackend reference");
    }
}

pub struct Mesh<B: MeshBackend> {
    pub backend: B,
}

impl<B: MeshBackend> Mesh<B> {
    pub fn num_vertices(&self) -> usize { self.backend.num_vertices() }
    pub fn num_edges(&self) -> usize { self.backend.num_edges() }
    pub fn num_faces(&self) -> usize { self.backend.num_faces() }
    pub fn num_halfedges(&self) -> usize { self.backend.num_halfedges() }
    pub fn num_corners(&self) -> usize { self.backend.num_corners() }
    pub fn num_boundaries(&self) -> usize { self.backend.num_boundaries() }

    pub fn euler_characteristic(&self) -> i32 {
        (self.num_vertices() as i32) - (self.num_edges() as i32) + (self.num_faces() as i32)
    }

    pub fn has_isolated_vertices(&self) -> bool {
        for v_idx in 0..self.num_vertices() {
            if self.backend.vertex_halfedge(v_idx).is_none() {
                eprintln!("Mesh has isolated vertices!");
                return true;
            }
        }
        false
    }

    pub fn has_isolated_faces(&self) -> bool {
        for f_idx in 0..self.num_faces() {
            let mut boundary_edges = 0;
            for h_idx in self.face_adjacent_halfedges(f_idx, true) {
                let twin_idx = self.backend.halfedge_twin(h_idx).expect("Twin should exist");
                if self.backend.halfedge_on_boundary(twin_idx) {
                    boundary_edges += 1;
                }
            }
            if boundary_edges == 3 {
                eprintln!("Mesh has isolated faces!");
                return true;
            }
        }
        false
    }

    pub fn has_non_manifold_vertices(&self) -> bool {
        let n_v = self.num_vertices();
        let mut adjacent_faces = vec![0; n_v];
        for f_idx in 0..self.num_faces() {
            for v_idx in self.face_adjacent_vertices(f_idx, true) {
                adjacent_faces[v_idx] += 1;
            }
        }
        for b_idx in 0..self.num_boundaries() {
            for v_idx in self.boundary_adjacent_vertices(b_idx, true) {
                adjacent_faces[v_idx] += 1;
            }
        }
        for (v_idx, &count) in adjacent_faces.iter().enumerate() {
            if count != self.vertex_degree(v_idx) {
                eprintln!("Mesh has non-manifold vertices!");
                return true;
            }
        }
        false
    }

    pub fn vertex_degree(&self, v_idx: usize) -> usize {
        self.vertex_adjacent_edges(v_idx, true).count()
    }

    pub fn is_isolated(&self, v_idx: usize) -> bool {
        self.backend.vertex_halfedge(v_idx).is_none()
    }

    pub fn on_boundary(&self, v_idx: usize) -> bool {
        for h_idx in self.vertex_adjacent_halfedges(v_idx, true) {
            if self.backend.halfedge_on_boundary(h_idx) {
                return true;
            }
        }
        false
    }

    // Iterators
    pub fn vertex_adjacent_halfedges(&self, v_idx: usize, ccw: bool) -> Box<dyn Iterator<Item = usize> + '_> {
        let start = self.backend.vertex_halfedge(v_idx);
        let mut curr = start;
        let mut just_started = true;
        Box::new(std::iter::from_fn(move || {
            if let Some(c) = curr {
                if !just_started && Some(c) == start { return None; }
                just_started = false;
                let val = c;
                if ccw {
                    if let Some(twin_idx) = self.backend.halfedge_twin(c) {
                        curr = self.backend.halfedge_next(twin_idx);
                    } else {
                        curr = None; // Boundary reached
                    }
                } else {
                    if let Some(prev_idx) = self.backend.halfedge_prev(c) {
                        curr = self.backend.halfedge_twin(prev_idx);
                    } else {
                        curr = None;
                    }
                }
                Some(val)
            } else { None }
        }))
    }

    pub fn vertex_adjacent_vertices(&self, v_idx: usize, ccw: bool) -> Box<dyn Iterator<Item = usize> + '_> {
        Box::new(self.vertex_adjacent_halfedges(v_idx, ccw).map(move |h| {
            let next_idx = self.backend.halfedge_next(h).expect("Next should exist");
            self.backend.halfedge_vertex(next_idx).expect("Vertex should exist")
        }))
    }

    pub fn vertex_adjacent_edges(&self, v_idx: usize, ccw: bool) -> Box<dyn Iterator<Item = usize> + '_> {
        Box::new(self.vertex_adjacent_halfedges(v_idx, ccw).map(move |h| self.backend.halfedge_edge(h).expect("Edge should exist")))
    }

    pub fn vertex_adjacent_faces(&self, v_idx: usize, ccw: bool) -> Box<dyn Iterator<Item = usize> + '_> {
        Box::new(self.vertex_adjacent_halfedges(v_idx, ccw).filter_map(move |h| {
            if !self.backend.halfedge_on_boundary(h) { self.backend.halfedge_face(h) } else { None }
        }))
    }

    pub fn vertex_adjacent_corners(&self, v_idx: usize, ccw: bool) -> Box<dyn Iterator<Item = usize> + '_> {
        Box::new(self.vertex_adjacent_halfedges(v_idx, ccw).filter_map(move |h_idx| {
            if !self.backend.halfedge_on_boundary(h_idx) {
                let next_idx = self.backend.halfedge_next(h_idx).expect("Next should exist");
                self.backend.halfedge_corner(next_idx)
            } else {
                None
            }
        }))
    }

    pub fn face_adjacent_halfedges(&self, f_idx: usize, ccw: bool) -> Box<dyn Iterator<Item = usize> + '_> {
        let start = self.backend.face_halfedge(f_idx);
        let mut curr = start;
        let mut just_started = true;
        Box::new(std::iter::from_fn(move || {
            if let Some(c) = curr {
                if !just_started && Some(c) == start { return None; }
                just_started = false;
                let val = c;
                curr = if ccw { self.backend.halfedge_next(c) } else { self.backend.halfedge_prev(c) };
                Some(val)
            } else { None }
        }))
    }

    pub fn face_adjacent_vertices(&self, f_idx: usize, ccw: bool) -> Box<dyn Iterator<Item = usize> + '_> {
        Box::new(self.face_adjacent_halfedges(f_idx, ccw).map(move |h| self.backend.halfedge_vertex(h).expect("Vertex should exist")))
    }

    pub fn face_adjacent_edges(&self, f_idx: usize, ccw: bool) -> Box<dyn Iterator<Item = usize> + '_> {
        Box::new(self.face_adjacent_halfedges(f_idx, ccw).map(move |h| self.backend.halfedge_edge(h).expect("Edge should exist")))
    }

    pub fn face_adjacent_faces(&self, f_idx: usize, ccw: bool) -> Box<dyn Iterator<Item = usize> + '_> {
        Box::new(self.face_adjacent_halfedges(f_idx, ccw).filter_map(move |h| {
            let twin_idx = self.backend.halfedge_twin(h).expect("Twin should exist");
            if !self.backend.halfedge_on_boundary(twin_idx) { self.backend.halfedge_face(twin_idx) } else { None }
        }))
    }

    pub fn face_adjacent_corners(&self, f_idx: usize, ccw: bool) -> Box<dyn Iterator<Item = usize> + '_> {
        Box::new(self.face_adjacent_halfedges(f_idx, ccw).filter_map(move |h| self.backend.halfedge_corner(h)))
    }

    pub fn face_adjacent_vertices_boundary(&self, f_idx: usize, ccw: bool) -> Box<dyn Iterator<Item = usize> + '_> {
        Box::new(self.face_adjacent_halfedges(f_idx, ccw).map(move |h| self.backend.halfedge_vertex(h).expect("Vertex should exist")))
    }

    pub fn boundary_adjacent_halfedges(&self, b_idx: usize, ccw: bool) -> Box<dyn Iterator<Item = usize> + '_> {
        let start = self.backend.boundary_halfedge(b_idx);
        let mut curr = start;
        let mut just_started = true;
        Box::new(std::iter::from_fn(move || {
            if let Some(c) = curr {
                if !just_started && Some(c) == start { return None; }
                just_started = false;
                let val = c;
                curr = if ccw { self.backend.halfedge_next(c) } else { self.backend.halfedge_prev(c) };
                Some(val)
            } else { None }
        }))
    }

    pub fn boundary_adjacent_vertices(&self, b_idx: usize, ccw: bool) -> Box<dyn Iterator<Item = usize> + '_> {
        Box::new(self.boundary_adjacent_halfedges(b_idx, ccw).map(move |h| self.backend.halfedge_vertex(h).expect("Vertex should exist")))
    }

    pub fn boundary_adjacent_edges(&self, b_idx: usize, ccw: bool) -> Box<dyn Iterator<Item = usize> + '_> {
        Box::new(self.boundary_adjacent_halfedges(b_idx, ccw).map(move |h| self.backend.halfedge_edge(h).expect("Edge should exist")))
    }
}

pub struct HalfEdgeBackend {
    pub vertex_halfedge: Vec<Option<usize>>,
    pub edge_halfedge: Vec<Option<usize>>,
    pub face_halfedge: Vec<Option<usize>>,
    pub corner_halfedge: Vec<Option<usize>>,
    pub boundary_halfedge: Vec<Option<usize>>,

    pub halfedge_vertex: Vec<Option<usize>>,
    pub halfedge_edge: Vec<Option<usize>>,
    pub halfedge_face: Vec<Option<usize>>,
    pub halfedge_corner: Vec<Option<usize>>,
    pub halfedge_next: Vec<Option<usize>>,
    pub halfedge_prev: Vec<Option<usize>>,
    pub halfedge_twin: Vec<Option<usize>>,
    pub halfedge_on_boundary: Vec<bool>,

    pub generators: Vec<Vec<usize>>,
}

impl MeshBackend for HalfEdgeBackend {
    fn num_vertices(&self) -> usize { self.vertex_halfedge.len() }
    fn num_edges(&self) -> usize { self.edge_halfedge.len() }
    fn num_faces(&self) -> usize { self.face_halfedge.len() }
    fn num_halfedges(&self) -> usize { self.halfedge_vertex.len() }
    fn num_corners(&self) -> usize { self.corner_halfedge.len() }
    fn num_boundaries(&self) -> usize { self.boundary_halfedge.len() }
    fn num_generators(&self) -> usize { self.generators.len() }

    fn vertex_halfedge(&self, v_idx: usize) -> Option<usize> { self.vertex_halfedge[v_idx] }
    fn edge_halfedge(&self, e_idx: usize) -> Option<usize> { self.edge_halfedge[e_idx] }
    fn face_halfedge(&self, f_idx: usize) -> Option<usize> { self.face_halfedge[f_idx] }
    fn corner_halfedge(&self, c_idx: usize) -> Option<usize> { self.corner_halfedge[c_idx] }
    fn boundary_halfedge(&self, b_idx: usize) -> Option<usize> { self.boundary_halfedge[b_idx] }
    fn generator(&self, g_idx: usize) -> &[usize] { &self.generators[g_idx] }

    fn halfedge_vertex(&self, h_idx: usize) -> Option<usize> { self.halfedge_vertex[h_idx] }
    fn halfedge_edge(&self, h_idx: usize) -> Option<usize> { self.halfedge_edge[h_idx] }
    fn halfedge_face(&self, h_idx: usize) -> Option<usize> { self.halfedge_face[h_idx] }
    fn halfedge_corner(&self, h_idx: usize) -> Option<usize> { self.halfedge_corner[h_idx] }
    fn halfedge_next(&self, h_idx: usize) -> Option<usize> { self.halfedge_next[h_idx] }
    fn halfedge_prev(&self, h_idx: usize) -> Option<usize> { self.halfedge_prev[h_idx] }
    fn halfedge_twin(&self, h_idx: usize) -> Option<usize> { self.halfedge_twin[h_idx] }
    fn halfedge_on_boundary(&self, h_idx: usize) -> bool { self.halfedge_on_boundary[h_idx] }

    fn add_generator(&mut self, generator: Vec<usize>) {
        self.generators.push(generator);
    }
}

#[derive(Clone)]
pub struct PolygonSoup {
    pub v: Vec<Vector>,
    pub f: Vec<usize>,
}

impl HalfEdgeBackend {
    pub fn new() -> Self {
        Self {
            vertex_halfedge: Vec::new(),
            edge_halfedge: Vec::new(),
            face_halfedge: Vec::new(),
            corner_halfedge: Vec::new(),
            boundary_halfedge: Vec::new(),
            halfedge_vertex: Vec::new(),
            halfedge_edge: Vec::new(),
            halfedge_face: Vec::new(),
            halfedge_corner: Vec::new(),
            halfedge_next: Vec::new(),
            halfedge_prev: Vec::new(),
            halfedge_twin: Vec::new(),
            halfedge_on_boundary: Vec::new(),
            generators: Vec::new(),
        }
    }
}

pub type DefaultMesh = Mesh<HalfEdgeBackend>;

impl Mesh<HalfEdgeBackend> {
    pub fn new() -> Self {
        Self {
            backend: HalfEdgeBackend::new(),
        }
    }

    pub fn build(&mut self, polygon_soup: &PolygonSoup) -> bool {
        let positions = &polygon_soup.v;
        let indices = &polygon_soup.f;
        
        let mut n_boundary_halfedges = 0;
        let mut sorted_edges = HashMap::new();
        for chunk in indices.chunks(3) {
            for j in 0..3 {
                let k = (j + 1) % 3;
                let mut u = chunk[j];
                let mut v = chunk[k];
                if u > v { std::mem::swap(&mut u, &mut v); }
                let key = (u, v);
                if sorted_edges.contains_key(&key) {
                    n_boundary_halfedges -= 1;
                } else {
                    sorted_edges.insert(key, true);
                    n_boundary_halfedges += 1;
                }
            }
        }

        let n_edges = sorted_edges.len();
        let n_faces = indices.len() / 3;
        let n_halfedges = 2 * n_edges;
        let n_interior_halfedges = n_halfedges - n_boundary_halfedges;

        // Idiomatic preallocation using vec!
        self.backend.vertex_halfedge = vec![None; positions.len()];
        self.backend.edge_halfedge = vec![None; n_edges];
        self.backend.face_halfedge = vec![None; n_faces];
        self.backend.halfedge_vertex = vec![None; n_halfedges];
        self.backend.halfedge_edge = vec![None; n_halfedges];
        self.backend.halfedge_face = vec![None; n_halfedges];
        self.backend.halfedge_corner = vec![None; n_halfedges];
        self.backend.halfedge_next = vec![None; n_halfedges];
        self.backend.halfedge_prev = vec![None; n_halfedges];
        self.backend.halfedge_twin = vec![None; n_halfedges];
        self.backend.halfedge_on_boundary = vec![false; n_halfedges];
        self.backend.corner_halfedge = vec![None; n_interior_halfedges];
        self.backend.boundary_halfedge.clear();
        self.backend.generators.clear();

        let mut has_twin_halfedge = vec![false; n_halfedges];
        let mut edge_count = HashMap::new();
        let mut existing_halfedges = HashMap::new();
        let mut e_index = 0;

        for (f_idx, chunk) in indices.chunks(3).enumerate() {
            for j in 0..3 {
                let k = (j + 1) % 3;
                let i_vert = chunk[j];
                let j_vert = chunk[k];

                let h_idx = f_idx * 3 + j;
                let next_idx = f_idx * 3 + k;
                let prev_idx = f_idx * 3 + (j + 2) % 3;

                self.backend.halfedge_next[h_idx] = Some(next_idx);
                self.backend.halfedge_prev[h_idx] = Some(prev_idx);
                self.backend.halfedge_on_boundary[h_idx] = false;

                self.backend.halfedge_vertex[h_idx] = Some(i_vert);
                self.backend.vertex_halfedge[i_vert] = Some(h_idx);

                self.backend.halfedge_face[h_idx] = Some(f_idx);
                self.backend.face_halfedge[f_idx] = Some(h_idx);

                let (mut u, mut v) = (i_vert, j_vert);
                if u > v { std::mem::swap(&mut u, &mut v); }
                let key = (u, v);

                if let Some(&twin_idx) = existing_halfedges.get(&key) {
                    self.backend.halfedge_twin[h_idx] = Some(twin_idx);
                    self.backend.halfedge_twin[twin_idx] = Some(h_idx);
                    let e_idx = self.backend.halfedge_edge[twin_idx].expect("Edge should exist");
                    self.backend.halfedge_edge[h_idx] = Some(e_idx);

                    has_twin_halfedge[h_idx] = true;
                    has_twin_halfedge[twin_idx] = true;
                    *edge_count.entry(key).or_insert(0) += 1;
                } else {
                    let e_idx = e_index;
                    e_index += 1;
                    self.backend.edge_halfedge[e_idx] = Some(h_idx);
                    self.backend.halfedge_edge[h_idx] = Some(e_idx);

                    existing_halfedges.insert(key, h_idx);
                    edge_count.insert(key, 1);
                }

                if *edge_count.get(&key).unwrap() > 2 {
                    eprintln!("Mesh has non-manifold edges!");
                    return false;
                }
            }
        }

        let mut h_idx_boundary = indices.len();
        let mut c_index = 0;
        for i in 0..indices.len() {
            let h_idx = i;
            if !has_twin_halfedge[h_idx] {
                let _b_idx = self.backend.boundary_halfedge.len();
                
                let mut boundary_cycle = Vec::new();
                let mut curr_h = h_idx;
                loop {
                    let bh_idx = h_idx_boundary;
                    h_idx_boundary += 1;
                    boundary_cycle.push(bh_idx);

                    let mut next_he = self.backend.halfedge_next[curr_h].expect("Next should exist");
                    while has_twin_halfedge[next_he] {
                        let twin_idx = self.backend.halfedge_twin[next_he].expect("Twin should exist");
                        next_he = self.backend.halfedge_next[twin_idx].expect("Next should exist");
                    }

                    self.backend.halfedge_vertex[bh_idx] = self.backend.halfedge_vertex[next_he];
                    self.backend.halfedge_edge[bh_idx] = self.backend.halfedge_edge[curr_h];
                    self.backend.halfedge_on_boundary[bh_idx] = true;
                    self.backend.halfedge_face[bh_idx] = None; // Boundaries usually don't have a face index in some parts, or use separate index.
                    // But our MeshBackend has boundary_halfedge.
                    
                    self.backend.halfedge_twin[bh_idx] = Some(curr_h);
                    self.backend.halfedge_twin[curr_h] = Some(bh_idx);

                    curr_h = next_he;
                    if curr_h == h_idx { break; }
                }

                let n = boundary_cycle.len();
                for j in 0..n {
                    let bh_curr = boundary_cycle[j];
                    self.backend.halfedge_next[bh_curr] = Some(boundary_cycle[(j + n - 1) % n]);
                    self.backend.halfedge_prev[bh_curr] = Some(boundary_cycle[(j + 1) % n]);
                    
                    has_twin_halfedge[bh_curr] = true;
                    has_twin_halfedge[self.backend.halfedge_twin[bh_curr].unwrap()] = true;
                }
                self.backend.boundary_halfedge.push(Some(boundary_cycle[0]));
            }

            if !self.backend.halfedge_on_boundary[h_idx] {
                self.backend.corner_halfedge[c_index] = Some(h_idx);
                self.backend.halfedge_corner[h_idx] = Some(c_index);
                c_index += 1;
            }
        }

        if self.has_isolated_vertices() || self.has_isolated_faces() || self.has_non_manifold_vertices() {
            return false;
        }

        true
    }
}
