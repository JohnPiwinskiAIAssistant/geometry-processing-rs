use crate::core::vertex::Vertex;
use crate::core::edge::Edge;
use crate::core::face::Face;
use crate::core::halfedge::Halfedge;
use crate::core::corner::Corner;
use crate::linear_algebra::Vector;
use std::collections::HashMap;

pub struct Mesh {
    pub vertices: Vec<Vertex>,
    pub edges: Vec<Edge>,
    pub faces: Vec<Face>,
    pub corners: Vec<Corner>,
    pub halfedges: Vec<Halfedge>,
    pub boundaries: Vec<Face>,
    pub generators: Vec<Vec<usize>>,
}

#[derive(Clone)]
pub struct PolygonSoup {
    pub v: Vec<Vector>,
    pub f: Vec<usize>,
}

impl Mesh {
    pub fn new() -> Self {
        Self {
            vertices: Vec::new(),
            edges: Vec::new(),
            faces: Vec::new(),
            corners: Vec::new(),
            halfedges: Vec::new(),
            boundaries: Vec::new(),
            generators: Vec::new(),
        }
    }

    pub fn euler_characteristic(&self) -> i32 {
        (self.vertices.len() as i32) - (self.edges.len() as i32) + (self.faces.len() as i32)
    }

    pub fn build(&mut self, polygon_soup: &PolygonSoup) -> bool {
        let positions = &polygon_soup.v;
        let indices = &polygon_soup.f;
        
        self.preallocate_elements(positions.len(), indices);

        let mut has_twin_halfedge = vec![false; self.halfedges.len()];
        let mut edge_count = HashMap::new();
        let mut existing_halfedges = HashMap::new();
        let mut e_index = 0;

        for (f_idx, chunk) in indices.chunks(3).enumerate() {
            self.faces[f_idx] = Face::new(f_idx);

            for j in 0..3 {
                let h_idx = f_idx * 3 + j;
                self.halfedges[h_idx] = Halfedge::new(h_idx);
            }

            for j in 0..3 {
                let k = (j + 1) % 3;
                let i_vert = chunk[j];
                let j_vert = chunk[k];

                let h_idx = f_idx * 3 + j;
                let next_idx = f_idx * 3 + k;
                let prev_idx = f_idx * 3 + (j + 2) % 3;

                self.halfedges[h_idx].next = Some(next_idx);
                self.halfedges[h_idx].prev = Some(prev_idx);
                self.halfedges[h_idx].on_boundary = false;

                self.halfedges[h_idx].vertex = Some(i_vert);
                self.vertices[i_vert].halfedge = Some(h_idx);

                self.halfedges[h_idx].face = Some(f_idx);
                self.faces[f_idx].halfedge = Some(h_idx);

                let (mut u, mut v) = (i_vert, j_vert);
                if u > v { std::mem::swap(&mut u, &mut v); }
                let key = (u, v);

                if let Some(&twin_idx) = existing_halfedges.get(&key) {
                    self.halfedges[h_idx].twin = Some(twin_idx);
                    self.halfedges[twin_idx].twin = Some(h_idx);
                    let e_idx = self.halfedges[twin_idx].edge.expect("Edge should exist");
                    self.halfedges[h_idx].edge = Some(e_idx);

                    has_twin_halfedge[h_idx] = true;
                    has_twin_halfedge[twin_idx] = true;
                    *edge_count.entry(key).or_insert(0) += 1;
                } else {
                    let e_idx = e_index;
                    e_index += 1;
                    self.edges[e_idx] = Edge::new(e_idx);
                    self.halfedges[h_idx].edge = Some(e_idx);
                    self.edges[e_idx].halfedge = Some(h_idx);

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
                let f_idx = self.boundaries.len();
                let mut b_face = Face::new(f_idx);

                let mut boundary_cycle = Vec::new();
                let mut curr_h = h_idx;
                loop {
                    let bh_idx = h_idx_boundary;
                    h_idx_boundary += 1;
                    let mut bh = Halfedge::new(bh_idx);
                    boundary_cycle.push(bh_idx);

                    let mut next_he = self.halfedges[curr_h].next.expect("Next should exist");
                    while has_twin_halfedge[next_he] {
                        let twin_idx = self.halfedges[next_he].twin.expect("Twin should exist");
                        next_he = self.halfedges[twin_idx].next.expect("Next should exist");
                    }

                    bh.vertex = self.halfedges[next_he].vertex;
                    bh.edge = self.halfedges[curr_h].edge;
                    bh.on_boundary = true;
                    bh.face = Some(f_idx);
                    
                    self.halfedges[bh_idx] = bh;
                    self.halfedges[bh_idx].twin = Some(curr_h);
                    self.halfedges[curr_h].twin = Some(bh_idx);

                    curr_h = next_he;
                    if curr_h == h_idx { break; }
                }

                let n = boundary_cycle.len();
                for j in 0..n {
                    let bh_curr = boundary_cycle[j];
                    self.halfedges[bh_curr].next = Some(boundary_cycle[(j + n - 1) % n]);
                    self.halfedges[bh_curr].prev = Some(boundary_cycle[(j + 1) % n]);
                    
                    has_twin_halfedge[bh_curr] = true;
                    has_twin_halfedge[self.halfedges[bh_curr].twin.unwrap()] = true;
                }
                b_face.halfedge = Some(boundary_cycle[0]);
                self.boundaries.push(b_face);
            }

            if !self.halfedges[h_idx].on_boundary {
                self.corners[c_index].halfedge = Some(h_idx);
                self.halfedges[h_idx].corner = Some(c_index);
                c_index += 1;
            }
        }

        if self.has_isolated_vertices() || self.has_isolated_faces() || self.has_non_manifold_vertices() {
            return false;
        }

        self.index_elements();
        true
    }

    fn preallocate_elements(&mut self, n_positions: usize, indices: &[usize]) {
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

        self.vertices = (0..n_positions).map(Vertex::new).collect();
        self.edges = (0..n_edges).map(Edge::new).collect();
        self.faces = (0..n_faces).map(Face::new).collect();
        self.halfedges = (0..n_halfedges).map(Halfedge::new).collect();
        self.corners = (0..n_interior_halfedges).map(Corner::new).collect();
        self.boundaries.clear();
        self.generators.clear();
    }

    pub fn has_isolated_vertices(&self) -> bool {
        for v in &self.vertices {
            if v.halfedge.is_none() {
                eprintln!("Mesh has isolated vertices!");
                return true;
            }
        }
        false
    }

    pub fn has_isolated_faces(&self) -> bool {
        for f in &self.faces {
            let mut boundary_edges = 0;
            for h_idx in self.face_adjacent_halfedges(f.index, true) {
                let twin_idx = self.halfedges[h_idx].twin.expect("Twin should exist");
                if self.halfedges[twin_idx].on_boundary {
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
        let mut adjacent_faces = vec![0; self.vertices.len()];
        for f in &self.faces {
            for v_idx in self.face_adjacent_vertices(f.index, true) {
                adjacent_faces[v_idx] += 1;
            }
        }
        for (b_idx, _) in self.boundaries.iter().enumerate() {
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

    pub fn index_elements(&mut self) {
        for (i, v) in self.vertices.iter_mut().enumerate() { v.index = i; }
        for (i, e) in self.edges.iter_mut().enumerate() { e.index = i; }
        for (i, f) in self.faces.iter_mut().enumerate() { f.index = i; }
        for (i, h) in self.halfedges.iter_mut().enumerate() { h.index = i; }
        for (i, c) in self.corners.iter_mut().enumerate() { c.index = i; }
        for (i, b) in self.boundaries.iter_mut().enumerate() { b.index = i; }
    }

    pub fn is_isolated(&self, v_idx: usize) -> bool {
        self.vertices[v_idx].halfedge.is_none()
    }

    pub fn on_boundary(&self, v_idx: usize) -> bool {
        for h_idx in self.vertex_adjacent_halfedges(v_idx, true) {
            if self.halfedges[h_idx].on_boundary {
                return true;
            }
        }
        false
    }

    // Iterators
    pub fn vertex_adjacent_halfedges(&self, v_idx: usize, ccw: bool) -> impl Iterator<Item = usize> + '_ {
        let start = self.vertices[v_idx].halfedge;
        let mut curr = start;
        let mut just_started = true;
        std::iter::from_fn(move || {
            if let Some(c) = curr {
                if !just_started && Some(c) == start { return None; }
                just_started = false;
                let val = c;
                let twin_idx = self.halfedges[c].twin.expect("Twin should exist");
                curr = if ccw {
                    self.halfedges[twin_idx].next
                } else {
                    let prev_idx = self.halfedges[c].prev.expect("Prev should exist");
                    self.halfedges[prev_idx].twin
                };
                Some(val)
            } else { None }
        })
    }

    pub fn vertex_adjacent_vertices(&self, v_idx: usize, ccw: bool) -> impl Iterator<Item = usize> + '_ {
        self.vertex_adjacent_halfedges(v_idx, ccw).map(move |h| {
            let twin_idx = self.halfedges[h].twin.expect("Twin should exist");
            self.halfedges[twin_idx].vertex.expect("Vertex should exist")
        })
    }

    pub fn vertex_adjacent_edges(&self, v_idx: usize, ccw: bool) -> impl Iterator<Item = usize> + '_ {
        self.vertex_adjacent_halfedges(v_idx, ccw).map(move |h| self.halfedges[h].edge.expect("Edge should exist"))
    }

    pub fn vertex_adjacent_faces(&self, v_idx: usize, ccw: bool) -> impl Iterator<Item = usize> + '_ {
        self.vertex_adjacent_halfedges(v_idx, ccw).filter_map(move |h| {
            if !self.halfedges[h].on_boundary { self.halfedges[h].face } else { None }
        })
    }

    pub fn vertex_adjacent_corners(&self, v_idx: usize, ccw: bool) -> impl Iterator<Item = usize> + '_ {
        self.vertex_adjacent_halfedges(v_idx, ccw).filter_map(move |h_idx| {
            if !self.halfedges[h_idx].on_boundary {
                let next_idx = self.halfedges[h_idx].next.expect("Next should exist");
                self.halfedges[next_idx].corner
            } else {
                None
            }
        })
    }

    pub fn face_adjacent_halfedges(&self, f_idx: usize, ccw: bool) -> impl Iterator<Item = usize> + '_ {
        let start = self.faces[f_idx].halfedge;
        let mut curr = start;
        let mut just_started = true;
        std::iter::from_fn(move || {
            if let Some(c) = curr {
                if !just_started && Some(c) == start { return None; }
                just_started = false;
                let val = c;
                curr = if ccw { self.halfedges[c].next } else { self.halfedges[c].prev };
                Some(val)
            } else { None }
        })
    }

    pub fn face_adjacent_vertices(&self, f_idx: usize, ccw: bool) -> impl Iterator<Item = usize> + '_ {
        self.face_adjacent_halfedges(f_idx, ccw).map(move |h| self.halfedges[h].vertex.expect("Vertex should exist"))
    }

    pub fn face_adjacent_edges(&self, f_idx: usize, ccw: bool) -> impl Iterator<Item = usize> + '_ {
        self.face_adjacent_halfedges(f_idx, ccw).map(move |h| self.halfedges[h].edge.expect("Edge should exist"))
    }

    pub fn face_adjacent_faces(&self, f_idx: usize, ccw: bool) -> impl Iterator<Item = usize> + '_ {
        self.face_adjacent_halfedges(f_idx, ccw).filter_map(move |h| {
            let twin_idx = self.halfedges[h].twin.expect("Twin should exist");
            if !self.halfedges[twin_idx].on_boundary { self.halfedges[twin_idx].face } else { None }
        })
    }

    pub fn face_adjacent_corners(&self, f_idx: usize, ccw: bool) -> impl Iterator<Item = usize> + '_ {
        self.face_adjacent_halfedges(f_idx, ccw).filter_map(move |h| self.halfedges[h].corner)
    }

    pub fn face_adjacent_vertices_boundary(&self, f_idx: usize, ccw: bool) -> impl Iterator<Item = usize> + '_ {
        self.face_adjacent_halfedges(f_idx, ccw).map(move |h| self.halfedges[h].vertex.expect("Vertex should exist"))
    }

    pub fn boundary_adjacent_vertices(&self, b_idx: usize, ccw: bool) -> impl Iterator<Item = usize> + '_ {
        let start = self.boundaries[b_idx].halfedge;
        let mut curr = start;
        let mut just_started = true;
        std::iter::from_fn(move || {
            if let Some(c) = curr {
                if !just_started && Some(c) == start { return None; }
                just_started = false;
                let val = self.halfedges[c].vertex.expect("Vertex should exist");
                curr = if ccw { self.halfedges[c].next } else { self.halfedges[c].prev };
                Some(val)
            } else { None }
        })
    }
}
