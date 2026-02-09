use crate::core::mesh::Mesh;
use crate::core::face::Face;
use crate::core::vertex::Vertex;
use crate::core::corner::Corner;

use crate::linear_algebra::{Vector, Complex, SparseMatrix, Triplet, DenseMatrix};
use crate::linear_algebra::traits::{SparseOps, DenseMatrixOps, Vector3Ops, Scalar};

pub struct Geometry<'a> {
    pub mesh: &'a Mesh,
    pub positions: Vec<Vector>,
}

impl<'a> Geometry<'a> {
    pub fn new(mesh: &'a Mesh, positions: Vec<Vector>, normalize_positions: bool) -> Self {
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
        let n = self.mesh.vertices.len();
        if n == 0 { return; }
        
        let mut cm = faer::Mat::zeros(3, 1);
        for i in 0..n {
            cm += &self.positions[i];
        }
        cm *= 1.0 / (n as f64);

        let mut radius = -1.0f64;
        for i in 0..n {
            self.positions[i] -= &cm;
            radius = radius.max((self.positions[i].transpose() * &self.positions[i])[(0, 0)].sqrt());
        }

        if rescale && radius > 0.0 {
            for i in 0..n {
                self.positions[i] *= 1.0 / radius;
            }
        }
    }

    pub fn vector(&self, h_idx: usize) -> Vector {
        let h = &self.mesh.halfedges[h_idx];
        let a = &self.positions[h.vertex.unwrap()];
        let b = &self.positions[self.mesh.halfedges[h.next.unwrap()].vertex.unwrap()];
        b - a
    }

    pub fn length(&self, e_idx: usize) -> f64 {
        (self.vector(self.mesh.edges[e_idx].halfedge.unwrap()).transpose() * &self.vector(self.mesh.edges[e_idx].halfedge.unwrap()))[(0, 0)].sqrt()
    }

    pub fn midpoint(&self, e_idx: usize) -> Vector {
        let e = &self.mesh.edges[e_idx];
        let h = &self.mesh.halfedges[e.halfedge.unwrap()];
        let twin = &self.mesh.halfedges[h.twin.unwrap()];
        let a = &self.positions[h.vertex.unwrap()];
        let b = &self.positions[twin.vertex.unwrap()];
        (a + b) * 0.5
    }

    pub fn mean_edge_length(&self) -> f64 {
        if self.mesh.edges.is_empty() { return 0.0; }
        let sum: f64 = (0..self.mesh.edges.len()).map(|i| self.length(i)).sum();
        sum / self.mesh.edges.len() as f64
    }

    pub fn area(&self, f: &Face) -> f64 {
        if self.mesh.halfedges[f.halfedge.unwrap()].on_boundary {
            return 0.0;
        }
        let h_idx = f.halfedge.unwrap();
        let u = self.vector(h_idx);
        let v = -&self.vector(self.mesh.halfedges[h_idx].prev.unwrap());
        
        // cross product for Mat<f64> (3x1)
        let cp = faer::mat![
            [u[(1, 0)] * v[(2, 0)] - u[(2, 0)] * v[(1, 0)]],
            [u[(2, 0)] * v[(0, 0)] - u[(0, 0)] * v[(2, 0)]],
            [u[(0, 0)] * v[(1, 0)] - u[(1, 0)] * v[(0, 0)]]
        ];
        0.5 * (cp.transpose() * &cp)[(0, 0)].sqrt()
    }

    pub fn total_area(&self) -> f64 {
        self.mesh.faces.iter().map(|f| self.area(f)).sum()
    }

    pub fn face_normal(&self, f: &Face) -> Option<Vector> {
        if self.mesh.halfedges[f.halfedge.unwrap()].on_boundary {
            return None;
        }
        let h_idx = f.halfedge.unwrap();
        let u = self.vector(h_idx);
        let v = -self.vector(self.mesh.halfedges[h_idx].prev.unwrap());
        
        let cp = u.cross(&v);
        let norm = cp.norm();
        Some(cp * (1.0 / norm))
    }

    pub fn centroid(&self, f: &Face) -> Vector {
        let h_idx = f.halfedge.unwrap();
        let a = &self.positions[self.mesh.halfedges[h_idx].vertex.unwrap()];
        let b = &self.positions[self.mesh.halfedges[self.mesh.halfedges[h_idx].next.unwrap()].vertex.unwrap()];
        let c = &self.positions[self.mesh.halfedges[self.mesh.halfedges[h_idx].prev.unwrap()].vertex.unwrap()];

        if self.mesh.halfedges[h_idx].on_boundary {
            return (a + b) * 0.5;
        }
        (a + b + c) * (1.0 / 3.0)
    }

    pub fn circumcenter(&self, f: &Face) -> Vector {
        let h_idx = f.halfedge.unwrap();
        let a = &self.positions[self.mesh.halfedges[h_idx].vertex.unwrap()];
        let b = &self.positions[self.mesh.halfedges[self.mesh.halfedges[h_idx].next.unwrap()].vertex.unwrap()];
        let c = &self.positions[self.mesh.halfedges[self.mesh.halfedges[h_idx].prev.unwrap()].vertex.unwrap()];

        if self.mesh.halfedges[h_idx].on_boundary {
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
        let e1 = self.vector(f.halfedge.unwrap()).unit();
        let normal = self.face_normal(f).unwrap();
        let e2 = normal.cross(&e1);
        [e1, e2]
    }

    pub fn angle(&self, c: &Corner) -> f64 {
        let h_idx = c.halfedge.unwrap();
        let u_raw = self.vector(self.mesh.halfedges[h_idx].prev.unwrap());
        let u = &u_raw * (1.0 / (u_raw.transpose() * &u_raw)[(0, 0)].sqrt());
        let v_raw = -&self.vector(self.mesh.halfedges[h_idx].next.unwrap());
        let v = &v_raw * (1.0 / (v_raw.transpose() * &v_raw)[(0, 0)].sqrt());
        let dot = (u.transpose() * v)[(0, 0)];
        dot.clamp(-1.0, 1.0).acos()
    }

    pub fn cotan(&self, h_idx: usize) -> f64 {
        let h = &self.mesh.halfedges[h_idx];
        if h.on_boundary {
            return 0.0;
        }

        let u = self.vector(h.prev.unwrap());
        let v = -&self.vector(h.next.unwrap());
        let dot = (u.transpose() * &v)[(0, 0)];
        let cp = faer::mat![
            [u[(1, 0)] * v[(2, 0)] - u[(2, 0)] * v[(1, 0)]],
            [u[(2, 0)] * v[(0, 0)] - u[(0, 0)] * v[(2, 0)]],
            [u[(0, 0)] * v[(1, 0)] - u[(1, 0)] * v[(0, 0)]]
        ];
        dot / (cp.transpose() * &cp)[(0, 0)].sqrt()
    }

    pub fn dihedral_angle(&self, h_idx: usize) -> f64 {
        let h = &self.mesh.halfedges[h_idx];
        let twin = &self.mesh.halfedges[h.twin.unwrap()];
        if h.on_boundary || twin.on_boundary {
            return 0.0;
        }

        let n1 = self.face_normal(&self.mesh.faces[h.face.unwrap()]).unwrap();
        let n2 = self.face_normal(&self.mesh.faces[twin.face.unwrap()]).unwrap();
        let w_raw = self.vector(h_idx);
        let w = &w_raw * (1.0 / (w_raw.transpose() * &w_raw)[(0, 0)].sqrt());

        let cos_theta = (n1.transpose() * &n2)[(0, 0)];
        let cp = faer::mat![
            [n1[(1, 0)] * n2[(2, 0)] - n1[(2, 0)] * n2[(1, 0)]],
            [n1[(2, 0)] * n2[(0, 0)] - n1[(0, 0)] * n2[(2, 0)]],
            [n1[(0, 0)] * n2[(1, 0)] - n1[(1, 0)] * n2[(0, 0)]]
        ];
        let sin_theta = (cp.transpose() * w)[(0, 0)];

        sin_theta.atan2(cos_theta)
    }

    pub fn barycentric_dual_area(&self, v: &Vertex) -> f64 {
        let mut area = 0.0;
        for f_idx in self.mesh.vertex_adjacent_faces(v.index, true) {
            area += self.area(&self.mesh.faces[f_idx]) / 3.0;
        }
        area
    }

    pub fn circumcentric_dual_area(&self, v: &Vertex) -> f64 {
        let mut area = 0.0;
        for h_idx in self.mesh.vertex_adjacent_halfedges(v.index, true) {
            let u2 = { let u = self.vector(self.mesh.halfedges[h_idx].prev.unwrap()); (u.transpose() * &u)[(0, 0)] };
            let v2 = { let v = self.vector(h_idx); (v.transpose() * &v)[(0, 0)] };
            let cot_alpha = self.cotan(self.mesh.halfedges[h_idx].prev.unwrap());
            let cot_beta = self.cotan(h_idx);

            area += (u2 * cot_alpha + v2 * cot_beta) / 8.0;
        }
        area
    }

    pub fn vertex_normal_equally_weighted(&self, v: &Vertex) -> Vector {
        let mut n = faer::Mat::zeros(3, 1);
        for f_idx in self.mesh.vertex_adjacent_faces(v.index, true) {
            n += &self.face_normal(&self.mesh.faces[f_idx]).unwrap();
        }
        let norm = (n.transpose() * &n)[(0, 0)].sqrt();
        &n * (1.0 / norm)
    }

    pub fn vertex_normal_area_weighted(&self, v: &Vertex) -> Vector {
        let mut n = faer::Mat::zeros(3, 1);
        for f_idx in self.mesh.vertex_adjacent_faces(v.index, true) {
            let normal = self.face_normal(&self.mesh.faces[f_idx]).unwrap();
            let area = self.area(&self.mesh.faces[f_idx]);
            n += &normal * area;
        }
        let norm = (n.transpose() * &n)[(0, 0)].sqrt();
        &n * (1.0 / norm)
    }

    pub fn vertex_normal_angle_weighted(&self, v: &Vertex) -> Vector {
        let mut n = faer::Mat::zeros(3, 1);
        for c_idx in self.mesh.vertex_adjacent_corners(v.index, true) {
            let c = &self.mesh.corners[c_idx];
            let normal = self.face_normal(&self.mesh.faces[self.mesh.halfedges[c.halfedge.unwrap()].face.unwrap()]).unwrap();
            let angle = self.angle(c);
            n += &normal * angle;
        }
        let norm = (n.transpose() * &n)[(0, 0)].sqrt();
        &n * (1.0 / norm)
    }

    pub fn vertex_normal_gauss_curvature(&self, v: &Vertex) -> Vector {
        let mut n = faer::Mat::zeros(3, 1);
        for h_idx in self.mesh.vertex_adjacent_halfedges(v.index, true) {
            let weight = 0.5 * self.dihedral_angle(h_idx) / self.length(self.mesh.halfedges[h_idx].edge.unwrap());
            n -= &self.vector(h_idx) * weight;
        }
        let norm = (n.transpose() * &n)[(0, 0)].sqrt();
        &n * (1.0 / norm)
    }

    pub fn vertex_normal_mean_curvature(&self, v: &Vertex) -> Vector {
        let mut n = faer::Mat::zeros(3, 1);
        for h_idx in self.mesh.vertex_adjacent_halfedges(v.index, true) {
            let weight = 0.5 * (self.cotan(h_idx) + self.cotan(self.mesh.halfedges[h_idx].twin.unwrap()));
            n -= &self.vector(h_idx) * weight;
        }
        let norm = (n.transpose() * &n)[(0, 0)].sqrt();
        &n * (1.0 / norm)
    }

    pub fn vertex_normal_sphere_inscribed(&self, v: &Vertex) -> Vector {
        let mut n = faer::Mat::zeros(3, 1);
        for c_idx in self.mesh.vertex_adjacent_corners(v.index, true) {
            let c = &self.mesh.corners[c_idx];
            let u = self.vector(self.mesh.halfedges[c.halfedge.unwrap()].prev.unwrap());
            let v_vec = -&self.vector(self.mesh.halfedges[c.halfedge.unwrap()].next.unwrap());
            
            let cp = faer::mat![
                [u[(1, 0)] * v_vec[(2, 0)] - u[(2, 0)] * v_vec[(1, 0)]],
                [u[(2, 0)] * v_vec[(0, 0)] - u[(0, 0)] * v_vec[(2, 0)]],
                [u[(0, 0)] * v_vec[(1, 0)] - u[(1, 0)] * v_vec[(0, 0)]]
            ];
            let u2 = (u.transpose() * &u)[(0, 0)];
            let v2 = (v_vec.transpose() * &v_vec)[(0, 0)];
            n += &cp * (1.0 / (u2 * v2));
        }
        let norm = (n.transpose() * &n)[(0, 0)].sqrt();
        &n * (1.0 / norm)
    }

    pub fn angle_defect(&self, v: &Vertex) -> f64 {
        let mut angle_sum = 0.0;
        for c_idx in self.mesh.vertex_adjacent_corners(v.index, true) {
            angle_sum += self.angle(&self.mesh.corners[c_idx]);
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
            sum += 0.5 * self.length(self.mesh.halfedges[h_idx].edge.expect("Edge should exist")) * self.dihedral_angle(h_idx);
        }
        sum
    }

    pub fn total_angle_defect(&self) -> f64 {
        self.mesh.vertices.iter().map(|v| self.angle_defect(v)).sum()
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
        let n = self.mesh.vertices.len();
        let mut t = Triplet::new(n, n);
        for v in &self.mesh.vertices {
            let i = v.index;
            let mut sum = 0.0;
            for h_idx in self.mesh.vertex_adjacent_halfedges(v.index, true) {
                let j = self.mesh.halfedges[self.mesh.halfedges[h_idx].twin.expect("Twin should exist")].vertex.expect("Vertex should exist");
                let weight = (self.cotan(h_idx) + self.cotan(self.mesh.halfedges[h_idx].twin.expect("Twin should exist"))) / 2.0;
                sum += weight;
                t.add_entry(S::from_f64(-weight), i, j);
            }
            t.add_entry(S::from_f64(sum), i, i);
        }
        SparseOps::from_triplets(n, n, &t.data)
    }

    pub fn laplace_matrix(&self) -> SparseMatrix<f64> {
        self.build_laplace_matrix::<f64>()
    }

    pub fn mass_matrix(&self) -> SparseMatrix<f64> {
        let n = self.mesh.vertices.len();
        let mut t = Triplet::new(n, n);
        for v in &self.mesh.vertices {
            t.add_entry(self.barycentric_dual_area(v), v.index, v.index);
        }
        SparseOps::from_triplets(n, n, &t.data)
    }

    pub fn complex_laplace_matrix(&self) -> SparseMatrix<Complex> {
        self.build_laplace_matrix::<Complex>()
    }
}
