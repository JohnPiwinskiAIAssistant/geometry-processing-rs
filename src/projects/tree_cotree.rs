use crate::core::mesh::{Mesh, MeshBackend};
use std::collections::{HashMap, VecDeque};

pub struct TreeCotree<'a, B: MeshBackend> {
    pub mesh: &'a mut Mesh<B>,
    pub vertex_parent: HashMap<usize, usize>,
    pub face_parent: HashMap<usize, usize>,
}

impl<'a, B: MeshBackend> TreeCotree<'a, B> {
    pub fn new(mesh: &'a mut Mesh<B>) -> Self {
        Self {
            mesh,
            vertex_parent: HashMap::new(),
            face_parent: HashMap::new(),
        }
    }

    fn build_primal_spanning_tree(&mut self) {
        let v_count = self.mesh.num_vertices();
        if v_count == 0 { return; }

        for i in 0..v_count {
            self.vertex_parent.insert(i, i);
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
        let u = self.mesh.backend.halfedge_vertex(h_idx).expect("Halfedge should have a vertex");
        let twin_idx = self.mesh.backend.halfedge_twin(h_idx).expect("Halfedge should have a twin");
        let v = self.mesh.backend.halfedge_vertex(twin_idx).expect("Twin should have a vertex");

        self.vertex_parent.get(&u) == Some(&v) || self.vertex_parent.get(&v) == Some(&u)
    }

    fn build_dual_spanning_cotree(&mut self) {
        let num_faces = self.mesh.num_faces();
        if num_faces == 0 { return; }

        for i in 0..num_faces {
            self.face_parent.insert(i, i);
        }

        let root = 0;
        let mut queue = VecDeque::new();
        queue.push_back(root);

        while let Some(f_idx) = queue.pop_front() {
            for h_idx in self.mesh.face_adjacent_halfedges(f_idx, true) {
                if !self.in_primal_spanning_tree(h_idx) {
                    let twin_idx = self.mesh.backend.halfedge_twin(h_idx).expect("Halfedge should have a twin");
                    if let Some(g_idx) = self.mesh.backend.halfedge_face(twin_idx) {
                        if !self.mesh.backend.halfedge_on_boundary(twin_idx) {
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
        let twin_idx = self.mesh.backend.halfedge_twin(h_idx).expect("Halfedge should have a twin");
        if let (Some(f), Some(g)) = (self.mesh.backend.halfedge_face(h_idx), self.mesh.backend.halfedge_face(twin_idx)) {
            if !self.mesh.backend.halfedge_on_boundary(h_idx) && !self.mesh.backend.halfedge_on_boundary(twin_idx) {
                return self.face_parent.get(&f) == Some(&g) || self.face_parent.get(&g) == Some(&f);
            }
        }
        false
    }

    fn shared_halfedge(&self, f: usize, g: usize) -> usize {
        for h_idx in self.mesh.face_adjacent_halfedges(f, true) {
            let twin_idx = self.mesh.backend.halfedge_twin(h_idx).expect("Halfedge should have a twin");
            if self.mesh.backend.halfedge_face(twin_idx) == Some(g) {
                return h_idx;
            }
        }
        panic!("Shared halfedge not found");
    }

    pub fn build_generators(&mut self) {
        self.build_primal_spanning_tree();
        self.build_dual_spanning_cotree();

        let num_edges = self.mesh.num_edges();
        for e_idx in 0..num_edges {
            let h_idx = self.mesh.backend.edge_halfedge(e_idx).expect("Edge should have a halfedge");
            if !self.in_primal_spanning_tree(h_idx) && !self.in_dual_spanning_tree(h_idx) {
                // trace back to root
                let mut temp_gen1 = Vec::new();
                let mut f = self.mesh.backend.halfedge_face(h_idx).expect("Halfedge should have a face");
                while let Some(&p) = self.face_parent.get(&f) {
                    if p == f { break; }
                    temp_gen1.push(self.shared_halfedge(f, p));
                    f = p;
                }

                let mut temp_gen2 = Vec::new();
                let twin_idx = self.mesh.backend.halfedge_twin(h_idx).expect("Halfedge should have a twin");
                let mut f2 = self.mesh.backend.halfedge_face(twin_idx).expect("Twin should have a face");
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
                    let twin = self.mesh.backend.halfedge_twin(temp_gen1[i]).expect("Twin should exist");
                    generator.push(twin);
                }
                for i in (0..=(n as usize)).rev() {
                    generator.push(temp_gen2[i]);
                }

                self.mesh.backend.add_generator(generator);
            }
        }
    }
}
