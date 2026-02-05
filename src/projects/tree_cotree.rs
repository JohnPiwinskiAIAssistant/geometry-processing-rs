use crate::core::mesh::Mesh;
use std::collections::{HashMap, VecDeque};

pub struct TreeCotree<'a> {
    pub mesh: &'a mut Mesh,
    pub vertex_parent: HashMap<usize, usize>,
    pub face_parent: HashMap<usize, usize>,
}

impl<'a> TreeCotree<'a> {
    pub fn new(mesh: &'a mut Mesh) -> Self {
        Self {
            mesh,
            vertex_parent: HashMap::new(),
            face_parent: HashMap::new(),
        }
    }

    fn build_primal_spanning_tree(&mut self) {
        let v_count = self.mesh.vertices.len();
        if v_count == 0 { return; }

        for v in &self.mesh.vertices {
            self.vertex_parent.insert(v.index, v.index);
        }

        let root = 0;
        let mut queue = VecDeque::new();
        queue.push_back(root);

        while let Some(u) = queue.pop_front() {
            for v in self.mesh.vertex_adjacent_vertices(u, true) {
                if let Some(&p) = self.vertex_parent.get(&v) {
                    if p == v && v != root {
                        self.vertex_parent.insert(v, u);
                        queue.push_back(v);
                    }
                }
            }
        }
    }

    fn in_primal_spanning_tree(&self, h_idx: usize) -> bool {
        let u = self.mesh.halfedges[h_idx].vertex.expect("Halfedge should have a vertex");
        let twin_idx = self.mesh.halfedges[h_idx].twin.expect("Halfedge should have a twin");
        let v = self.mesh.halfedges[twin_idx].vertex.expect("Twin should have a vertex");

        self.vertex_parent.get(&u) == Some(&v) || self.vertex_parent.get(&v) == Some(&u)
    }

    fn build_dual_spanning_cotree(&mut self) {
        if self.mesh.faces.is_empty() { return; }

        for f in &self.mesh.faces {
            self.face_parent.insert(f.index, f.index);
        }

        let root = 0;
        let mut queue = VecDeque::new();
        queue.push_back(root);

        while let Some(f_idx) = queue.pop_front() {
            for h_idx in self.mesh.face_adjacent_halfedges(f_idx, true) {
                if !self.in_primal_spanning_tree(h_idx) {
                    let twin_idx = self.mesh.halfedges[h_idx].twin.expect("Halfedge should have a twin");
                    if let Some(g_idx) = self.mesh.halfedges[twin_idx].face {
                        if !self.mesh.halfedges[twin_idx].on_boundary {
                            if let Some(&p) = self.face_parent.get(&g_idx) {
                                if p == g_idx && g_idx != root {
                                    self.face_parent.insert(g_idx, f_idx);
                                    queue.push_back(g_idx);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    fn in_dual_spanning_tree(&self, h_idx: usize) -> bool {
        let twin_idx = self.mesh.halfedges[h_idx].twin.expect("Halfedge should have a twin");
        if let (Some(f), Some(g)) = (self.mesh.halfedges[h_idx].face, self.mesh.halfedges[twin_idx].face) {
            if !self.mesh.halfedges[h_idx].on_boundary && !self.mesh.halfedges[twin_idx].on_boundary {
                return self.face_parent.get(&f) == Some(&g) || self.face_parent.get(&g) == Some(&f);
            }
        }
        false
    }

    fn shared_halfedge(&self, f: usize, g: usize) -> usize {
        for h_idx in self.mesh.face_adjacent_halfedges(f, true) {
            let twin_idx = self.mesh.halfedges[h_idx].twin.expect("Halfedge should have a twin");
            if self.mesh.halfedges[twin_idx].face == Some(g) {
                return h_idx;
            }
        }
        panic!("Shared halfedge not found");
    }

    pub fn build_generators(&mut self) {
        self.build_primal_spanning_tree();
        self.build_dual_spanning_cotree();

        for e_idx in 0..self.mesh.edges.len() {
            let h_idx = self.mesh.edges[e_idx].halfedge.expect("Edge should have a halfedge");
            if !self.in_primal_spanning_tree(h_idx) && !self.in_dual_spanning_tree(h_idx) {
                // trace back to root
                let mut temp_gen1 = Vec::new();
                let mut f = self.mesh.halfedges[h_idx].face.expect("Halfedge should have a face");
                while let Some(&p) = self.face_parent.get(&f) {
                    if p == f { break; }
                    temp_gen1.push(self.shared_halfedge(f, p));
                    f = p;
                }

                let mut temp_gen2 = Vec::new();
                let twin_idx = self.mesh.halfedges[h_idx].twin.expect("Halfedge should have a twin");
                let mut f2 = self.mesh.halfedges[twin_idx].face.expect("Twin should have a face");
                while let Some(&p) = self.face_parent.get(&f2) {
                    if p == f2 { break; }
                    temp_gen2.push(self.shared_halfedge(f2, p));
                    f2 = p;
                }

                let mut m = (temp_gen1.len() as i32) - 1;
                let mut n = (temp_gen2.len() as i32) - 1;
                while m >= 0 && n >= 0 && temp_gen1[m as usize] == temp_gen2[n as usize] {
                    m -= 1;
                    n -= 1;
                }

                let mut generator = vec![h_idx];
                for i in 0..=(m as usize) {
                    let twin = self.mesh.halfedges[temp_gen1[i]].twin.expect("Twin should exist");
                    generator.push(twin);
                }
                for i in (0..=(n as usize)).rev() {
                    generator.push(temp_gen2[i]);
                }

                self.mesh.generators.push(generator);
            }
        }
    }
}
