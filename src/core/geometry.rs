use crate::core::mesh::{Mesh, MeshBackend};
use crate::core::face::Face;
use crate::core::vertex::Vertex;
use crate::core::corner::Corner;

use crate::linear_algebra::{Vector, Complex, SparseMatrix, Triplet};
use crate::linear_algebra::traits::{SparseOps, DenseMatrixOps, Vector3Ops, Scalar};

pub struct Geometry<'a, B: MeshBackend> {
    pub mesh: &'a Mesh<B>,
    pub positions: Vec<Vector>,
}

impl<'a, B: MeshBackend> Geometry<'a, B> {
    pub fn new(mesh: &'a Mesh<B>, positions: Vec<Vector>, normalize_positions: bool) -> Self {
        let mut geom = Self {
            mesh,
            positions,
        };

        if normalize_positions {
            geom.normalize(true);
        }

        geom
    }

    pub fn normalize(&mut self, rescale: bool) {
        let n = self.mesh.num_vertices();
        if n == 0 { return; }
        
        let mut cm = faer::Mat::zeros(3, 1);
        for i in 0..n {
            cm += &self.positions[i];
        }
        cm *= 1.0 / (n as f64);

        let mut radius = -1.0f64;
        for i in 0..n {
            self.positions[i] -= &cm;
            radius = radius.max(self.positions[i].norm());
        }

        if rescale && radius > 0.0 {
            for i in 0..n {
                self.positions[i] *= 1.0 / radius;
            }
        }
    }

    pub fn vector(&self, h_idx: usize) -> Vector {
        let v_idx = self.mesh.backend.halfedge_vertex(h_idx).unwrap();
        let a = &self.positions[v_idx];
        let next_h = self.mesh.backend.halfedge_next(h_idx).unwrap();
        let next_v = self.mesh.backend.halfedge_vertex(next_h).unwrap();
        let b = &self.positions[next_v];
        b - a
    }

    pub fn length(&self, e_idx: usize) -> f64 {
        self.vector(self.mesh.backend.edge_halfedge(e_idx).unwrap()).norm()
    }

    pub fn midpoint(&self, e_idx: usize) -> Vector {
        let h_idx = self.mesh.backend.edge_halfedge(e_idx).unwrap();
        let twin_idx = self.mesh.backend.halfedge_twin(h_idx).unwrap();
        let a_v = self.mesh.backend.halfedge_vertex(h_idx).unwrap();
        let b_v = self.mesh.backend.halfedge_vertex(twin_idx).unwrap();
        let a = &self.positions[a_v];
        let b = &self.positions[b_v];
        (a + b) * 0.5
    }

    pub fn mean_edge_length(&self) -> f64 {
        let num_edges = self.mesh.num_edges();
        if num_edges == 0 { return 0.0; }
        let sum: f64 = (0..num_edges).map(|i| self.length(i)).sum();
        sum / num_edges as f64
    }

    pub fn area(&self, f: &Face) -> f64 {
        let h_idx = self.mesh.backend.face_halfedge(f.index).unwrap();
        if self.mesh.backend.halfedge_on_boundary(h_idx) {
            return 0.0;
        }
        let u = self.vector(h_idx);
        let prev_h = self.mesh.backend.halfedge_prev(h_idx).unwrap();
        let v = -&self.vector(prev_h);
        
        let cp = u.cross(&v);
        0.5 * cp.norm()
    }

    pub fn total_area(&self) -> f64 {
        (0..self.mesh.num_faces()).map(|i| self.area(&Face::new(i))).sum()
    }

    pub fn face_normal(&self, f: &Face) -> Option<Vector> {
        let h_idx = self.mesh.backend.face_halfedge(f.index).unwrap();
        if self.mesh.backend.halfedge_on_boundary(h_idx) {
            return None;
        }
        let u = self.vector(h_idx);
        let prev_h = self.mesh.backend.halfedge_prev(h_idx).unwrap();
        let v = -self.vector(prev_h);
        
        let cp = u.cross(&v);
        let norm = cp.norm();
        Some(cp * (1.0 / norm))
    }

    pub fn centroid(&self, f: &Face) -> Vector {
        let h_idx = self.mesh.backend.face_halfedge(f.index).unwrap();
        let next_h = self.mesh.backend.halfedge_next(h_idx).unwrap();
        let prev_h = self.mesh.backend.halfedge_prev(h_idx).unwrap();
        
        let a = &self.positions[self.mesh.backend.halfedge_vertex(h_idx).unwrap()];
        let b = &self.positions[self.mesh.backend.halfedge_vertex(next_h).unwrap()];
        let c = &self.positions[self.mesh.backend.halfedge_vertex(prev_h).unwrap()];

        if self.mesh.backend.halfedge_on_boundary(h_idx) {
            return (a + b) * 0.5;
        }
        (a + b + c) * (1.0 / 3.0)
    }

    pub fn circumcenter(&self, f: &Face) -> Vector {
        let h_idx = self.mesh.backend.face_halfedge(f.index).unwrap();
        let next_h = self.mesh.backend.halfedge_next(h_idx).unwrap();
        let prev_h = self.mesh.backend.halfedge_prev(h_idx).unwrap();

        let a = &self.positions[self.mesh.backend.halfedge_vertex(h_idx).unwrap()];
        let b = &self.positions[self.mesh.backend.halfedge_vertex(next_h).unwrap()];
        let c = &self.positions[self.mesh.backend.halfedge_vertex(prev_h).unwrap()];

        if self.mesh.backend.halfedge_on_boundary(h_idx) {
            return (a + b) * 0.5;
        }

        let ac = c - a;
        let ab = b - a;
        
        let w = ab.cross(&ac);

        let ac2 = ac.norm_sq();
        let ab2 = ab.norm_sq();
        let w2 = w.norm_sq();

        let w_x_ab = w.cross(&ab);
        let ac_x_w = ac.cross(&w);

        let u = w_x_ab * ac2;
        let v = ac_x_w * ab2;
        let x = (u + v) * (0.5 / w2);

        x + a
    }

    pub fn orthonormal_bases(&self, f: &Face) -> [Vector; 2] {
        let h_idx = self.mesh.backend.face_halfedge(f.index).unwrap();
        let e1 = self.vector(h_idx).unit();
        let normal = self.face_normal(f).unwrap();
        let e2 = normal.cross(&e1);
        [e1, e2]
    }

    pub fn angle(&self, c: &Corner) -> f64 {
        let h_idx = self.mesh.backend.corner_halfedge(c.index).expect("Corner should have a halfedge");
        let prev_h = self.mesh.backend.halfedge_prev(h_idx).expect("Prev should exist");
        let next_h = self.mesh.backend.halfedge_next(h_idx).expect("Next should exist");
        
        let u_raw = self.vector(prev_h);
        let u = &u_raw * (1.0 / u_raw.norm());
        let v_raw = -&self.vector(next_h);
        let v = &v_raw * (1.0 / v_raw.norm());
        let dot = u.dot(&v);
        dot.clamp(-1.0, 1.0).acos()
    }

    pub fn cotan(&self, h_idx: usize) -> f64 {
        if self.mesh.backend.halfedge_on_boundary(h_idx) {
            return 0.0;
        }

        let prev_h = self.mesh.backend.halfedge_prev(h_idx).unwrap();
        let next_h = self.mesh.backend.halfedge_next(h_idx).unwrap();
        let u = self.vector(prev_h);
        let v = -&self.vector(next_h);
        let dot = (u.transpose() * &v)[(0, 0)];
        let cp = u.cross(&v);
        dot / cp.norm()
    }

    pub fn dihedral_angle(&self, h_idx: usize) -> f64 {
        let twin_idx = self.mesh.backend.halfedge_twin(h_idx).unwrap();
        if self.mesh.backend.halfedge_on_boundary(h_idx) || self.mesh.backend.halfedge_on_boundary(twin_idx) {
            return 0.0;
        }

        let f_idx = self.mesh.backend.halfedge_face(h_idx).expect("Face should exist");
        let twin_f_idx = self.mesh.backend.halfedge_face(twin_idx).expect("Face should exist");
        
        let n1 = self.face_normal(&Face::new(f_idx)).unwrap();
        let n2 = self.face_normal(&Face::new(twin_f_idx)).unwrap();
        let w_raw = self.vector(h_idx);
        let w = &w_raw * (1.0 / w_raw.norm());

        let cos_theta = n1.dot(&n2);
        let cp = n1.cross(&n2);
        let sin_theta = cp.dot(&w);

        sin_theta.atan2(cos_theta)
    }

    pub fn barycentric_dual_area(&self, v: &Vertex) -> f64 {
        let mut area = 0.0;
        for f_idx in self.mesh.vertex_adjacent_faces(v.index, true) {
            area += self.area(&Face::new(f_idx)) / 3.0;
        }
        area
    }

    pub fn circumcentric_dual_area(&self, v: &Vertex) -> f64 {
        let mut area = 0.0;
        for h_idx in self.mesh.vertex_adjacent_halfedges(v.index, true) {
            let prev_h = self.mesh.backend.halfedge_prev(h_idx).unwrap();
            let u2 = self.vector(prev_h).norm_sq();
            let v2 = self.vector(h_idx).norm_sq();
            let cot_alpha = self.cotan(prev_h);
            let cot_beta = self.cotan(h_idx);

            area += (u2 * cot_alpha + v2 * cot_beta) / 8.0;
        }
        area
    }

    pub fn vertex_normal_equally_weighted(&self, v: &Vertex) -> Vector {
        let mut n = faer::Mat::zeros(3, 1);
        for f_idx in self.mesh.vertex_adjacent_faces(v.index, true) {
            n += &self.face_normal(&Face::new(f_idx)).unwrap();
        }
        let norm = n.norm();
        &n * (1.0 / norm)
    }

    pub fn vertex_normal_area_weighted(&self, v: &Vertex) -> Vector {
        let mut n = faer::Mat::zeros(3, 1);
        for f_idx in self.mesh.vertex_adjacent_faces(v.index, true) {
            let normal = self.face_normal(&Face::new(f_idx)).unwrap();
            let area = self.area(&Face::new(f_idx));
            n += &normal * area;
        }
        let norm = n.norm();
        &n * (1.0 / norm)
    }

    pub fn vertex_normal_angle_weighted(&self, v: &Vertex) -> Vector {
        let mut n = faer::Mat::zeros(3, 1);
        for c_idx in self.mesh.vertex_adjacent_corners(v.index, true) {
            let h_idx = self.mesh.backend.corner_halfedge(c_idx).expect("Corner should have a halfedge");
            let f_idx = self.mesh.backend.halfedge_face(h_idx).expect("Face should exist");
            let normal = self.face_normal(&Face::new(f_idx)).unwrap();
            let angle = self.angle(&Corner::new(c_idx));
            n += &normal * angle;
        }
        let norm = n.norm();
        &n * (1.0 / norm)
    }

    pub fn vertex_normal_gauss_curvature(&self, v: &Vertex) -> Vector {
        let mut n = faer::Mat::zeros(3, 1);
        for h_idx in self.mesh.vertex_adjacent_halfedges(v.index, true) {
            let e_idx = self.mesh.backend.halfedge_edge(h_idx).unwrap();
            let weight = 0.5 * self.dihedral_angle(h_idx) / self.length(e_idx);
            n -= &self.vector(h_idx) * weight;
        }
        let norm = n.norm();
        &n * (1.0 / norm)
    }

    pub fn vertex_normal_mean_curvature(&self, v: &Vertex) -> Vector {
        let mut n = faer::Mat::zeros(3, 1);
        for h_idx in self.mesh.vertex_adjacent_halfedges(v.index, true) {
            let twin_idx = self.mesh.backend.halfedge_twin(h_idx).unwrap();
            let weight = 0.5 * (self.cotan(h_idx) + self.cotan(twin_idx));
            n -= &self.vector(h_idx) * weight;
        }
        let norm = n.norm();
        &n * (1.0 / norm)
    }

    pub fn vertex_normal_sphere_inscribed(&self, v: &Vertex) -> Vector {
        let mut n = faer::Mat::zeros(3, 1);
        for c_idx in self.mesh.vertex_adjacent_corners(v.index, true) {
            let h_idx = self.mesh.backend.corner_halfedge(c_idx).expect("Corner should have a halfedge");
            let prev_h = self.mesh.backend.halfedge_prev(h_idx).unwrap();
            let next_h = self.mesh.backend.halfedge_next(h_idx).unwrap();
            let u = self.vector(prev_h);
            let v_vec = -&self.vector(next_h);
            
            let cp = u.cross(&v_vec);
            let u2 = u.norm_sq();
            let v2 = v_vec.norm_sq();
            n += &cp * (1.0 / (u2 * v2));
        }
        let norm = n.norm();
        &n * (1.0 / norm)
    }

    pub fn angle_defect(&self, v: &Vertex) -> f64 {
        let mut angle_sum = 0.0;
        for c_idx in self.mesh.vertex_adjacent_corners(v.index, true) {
            angle_sum += self.angle(&Corner::new(c_idx));
        }
        if self.mesh.on_boundary(v.index) {
            std::f64::consts::PI - angle_sum
        } else {
            2.0 * std::f64::consts::PI - angle_sum
        }
    }

    pub fn scalar_gauss_curvature(&self, v: &Vertex) -> f64 {
        self.angle_defect(v)
    }

    pub fn scalar_mean_curvature(&self, v: &Vertex) -> f64 {
        let mut sum = 0.0;
        for h_idx in self.mesh.vertex_adjacent_halfedges(v.index, true) {
            let e_idx = self.mesh.backend.halfedge_edge(h_idx).expect("Edge should exist");
            sum += 0.5 * self.length(e_idx) * self.dihedral_angle(h_idx);
        }
        sum
    }

    pub fn total_angle_defect(&self) -> f64 {
        (0..self.mesh.num_vertices()).map(|i| self.angle_defect(&Vertex::new(i))).sum()
    }

    pub fn principal_curvatures(&self, v: &Vertex) -> [f64; 2] {
        let a = self.circumcentric_dual_area(v);
        let h = self.scalar_mean_curvature(v) / a;
        let k = self.angle_defect(v) / a;

        let mut discriminant = h * h - k;
        if discriminant > 0.0 {
            discriminant = discriminant.sqrt();
        } else {
            discriminant = 0.0;
        }

        let k1 = h - discriminant;
        let k2 = h + discriminant;

        [k1, k2]
    }

    pub fn build_laplace_matrix<S: Scalar>(&self) -> SparseMatrix<S> {
        let n = self.mesh.num_vertices();
        let mut t = Triplet::new(n, n);
        for i in 0..n {
            let mut sum = 0.0;
            for h_idx in self.mesh.vertex_adjacent_halfedges(i, true) {
                let next_idx = self.mesh.backend.halfedge_next(h_idx).expect("Next should exist");
                let j = self.mesh.backend.halfedge_vertex(next_idx).expect("Vertex should exist");
                
                let weight = self.cotan(h_idx);
                let twin_weight = if let Some(twin_idx) = self.mesh.backend.halfedge_twin(h_idx) {
                    self.cotan(twin_idx)
                } else {
                    0.0
                };
                
                let total_weight = (weight + twin_weight) / 2.0;
                sum += total_weight;
                t.add_entry(S::from_f64(-total_weight), i, j);
            }
            t.add_entry(S::from_f64(sum), i, i);
        }
        SparseOps::from_triplets(n, n, &t.data)
    }

    pub fn laplace_matrix(&self) -> SparseMatrix<f64> {
        self.build_laplace_matrix::<f64>()
    }

    pub fn mass_matrix(&self) -> SparseMatrix<f64> {
        let n = self.mesh.num_vertices();
        let mut t = Triplet::new(n, n);
        for i in 0..n {
            t.add_entry(self.barycentric_dual_area(&Vertex::new(i)), i, i);
        }
        SparseOps::from_triplets(n, n, &t.data)
    }

    pub fn complex_laplace_matrix(&self) -> SparseMatrix<Complex> {
        self.build_laplace_matrix::<Complex>()
    }
}
