use crate::core::mesh::Mesh;
use crate::core::face::Face;
use crate::core::vertex::Vertex;
use crate::core::corner::Corner;

use crate::linear_algebra::{Vector, Complex, SparseMatrix, Triplet, ComplexTriplet, ComplexSparseMatrix};

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
        
        let mut cm = Vector::new(0.0, 0.0, 0.0);
        for i in 0..n {
            cm.increment_by(self.positions[i]);
        }
        cm.divide_by(n as f64);

        let mut radius = -1.0f64;
        for i in 0..n {
            self.positions[i].decrement_by(cm);
            radius = radius.max(self.positions[i].norm());
        }

        if rescale && radius > 0.0 {
            for i in 0..n {
                self.positions[i].divide_by(radius);
            }
        }
    }

    pub fn vector(&self, h_idx: usize) -> Vector {
        let h = &self.mesh.halfedges[h_idx];
        let a = self.positions[h.vertex.unwrap()];
        let b = self.positions[self.mesh.halfedges[h.next.unwrap()].vertex.unwrap()];
        b.minus(a)
    }

    pub fn length(&self, e_idx: usize) -> f64 {
        self.vector(self.mesh.edges[e_idx].halfedge.unwrap()).norm()
    }

    pub fn midpoint(&self, e_idx: usize) -> Vector {
        let e = &self.mesh.edges[e_idx];
        let h = &self.mesh.halfedges[e.halfedge.unwrap()];
        let twin = &self.mesh.halfedges[h.twin.unwrap()];
        let a = self.positions[h.vertex.unwrap()];
        let b = self.positions[twin.vertex.unwrap()];
        a.plus(b).over(2.0)
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
        let v = self.vector(self.mesh.halfedges[h_idx].prev.unwrap()).negated();
        0.5 * u.cross(v).norm()
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
        let v = self.vector(self.mesh.halfedges[h_idx].prev.unwrap()).negated();
        Some(u.cross(v).unit())
    }

    pub fn centroid(&self, f: &Face) -> Vector {
        let h_idx = f.halfedge.unwrap();
        let a = self.positions[self.mesh.halfedges[h_idx].vertex.unwrap()];
        let b = self.positions[self.mesh.halfedges[self.mesh.halfedges[h_idx].next.unwrap()].vertex.unwrap()];
        let c = self.positions[self.mesh.halfedges[self.mesh.halfedges[h_idx].prev.unwrap()].vertex.unwrap()];

        if self.mesh.halfedges[h_idx].on_boundary {
            return a.plus(b).over(2.0);
        }
        a.plus(b).plus(c).over(3.0)
    }

    pub fn circumcenter(&self, f: &Face) -> Vector {
        let h_idx = f.halfedge.unwrap();
        let a = self.positions[self.mesh.halfedges[h_idx].vertex.unwrap()];
        let b = self.positions[self.mesh.halfedges[self.mesh.halfedges[h_idx].next.unwrap()].vertex.unwrap()];
        let c = self.positions[self.mesh.halfedges[self.mesh.halfedges[h_idx].prev.unwrap()].vertex.unwrap()];

        if self.mesh.halfedges[h_idx].on_boundary {
            return a.plus(b).over(2.0);
        }

        let ac = c.minus(a);
        let ab = b.minus(a);
        let w = ab.cross(ac);

        let u = (w.cross(ab)).times(ac.norm2());
        let v = (ac.cross(w)).times(ab.norm2());
        let x = (u.plus(v)).over(2.0 * w.norm2());

        x.plus(a)
    }

    pub fn orthonormal_bases(&self, f: &Face) -> [Vector; 2] {
        let e1 = self.vector(f.halfedge.unwrap()).unit();
        let normal = self.face_normal(f).unwrap();
        let e2 = normal.cross(e1);
        [e1, e2]
    }

    pub fn angle(&self, c: &Corner) -> f64 {
        let h_idx = c.halfedge.unwrap();
        let u = self.vector(self.mesh.halfedges[h_idx].prev.unwrap()).unit();
        let v = self.vector(self.mesh.halfedges[h_idx].next.unwrap()).negated().unit();
        u.dot(v).clamp(-1.0, 1.0).acos()
    }

    pub fn cotan(&self, h_idx: usize) -> f64 {
        let h = &self.mesh.halfedges[h_idx];
        if h.on_boundary {
            return 0.0;
        }

        let u = self.vector(h.prev.unwrap());
        let v = self.vector(h.next.unwrap()).negated();
        u.dot(v) / u.cross(v).norm()
    }

    pub fn dihedral_angle(&self, h_idx: usize) -> f64 {
        let h = &self.mesh.halfedges[h_idx];
        let twin = &self.mesh.halfedges[h.twin.unwrap()];
        if h.on_boundary || twin.on_boundary {
            return 0.0;
        }

        let n1 = self.face_normal(&self.mesh.faces[h.face.unwrap()]).unwrap();
        let n2 = self.face_normal(&self.mesh.faces[twin.face.unwrap()]).unwrap();
        let w = self.vector(h_idx).unit();

        let cos_theta = n1.dot(n2);
        let sin_theta = n1.cross(n2).dot(w);

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
            let u2 = self.vector(self.mesh.halfedges[h_idx].prev.unwrap()).norm2();
            let v2 = self.vector(h_idx).norm2();
            let cot_alpha = self.cotan(self.mesh.halfedges[h_idx].prev.unwrap());
            let cot_beta = self.cotan(h_idx);

            area += (u2 * cot_alpha + v2 * cot_beta) / 8.0;
        }
        area
    }

    pub fn vertex_normal_equally_weighted(&self, v: &Vertex) -> Vector {
        let mut n = Vector::new(0.0, 0.0, 0.0);
        for f_idx in self.mesh.vertex_adjacent_faces(v.index, true) {
            n.increment_by(self.face_normal(&self.mesh.faces[f_idx]).unwrap());
        }
        n.unit()
    }

    pub fn vertex_normal_area_weighted(&self, v: &Vertex) -> Vector {
        let mut n = Vector::new(0.0, 0.0, 0.0);
        for f_idx in self.mesh.vertex_adjacent_faces(v.index, true) {
            let normal = self.face_normal(&self.mesh.faces[f_idx]).unwrap();
            let area = self.area(&self.mesh.faces[f_idx]);
            n.increment_by(normal.times(area));
        }
        n.unit()
    }

    pub fn vertex_normal_angle_weighted(&self, v: &Vertex) -> Vector {
        let mut n = Vector::new(0.0, 0.0, 0.0);
        for c_idx in self.mesh.vertex_adjacent_corners(v.index, true) {
            let c = &self.mesh.corners[c_idx];
            let normal = self.face_normal(&self.mesh.faces[self.mesh.halfedges[c.halfedge.unwrap()].face.unwrap()]).unwrap();
            let angle = self.angle(c);
            n.increment_by(normal.times(angle));
        }
        n.unit()
    }

    pub fn vertex_normal_gauss_curvature(&self, v: &Vertex) -> Vector {
        let mut n = Vector::new(0.0, 0.0, 0.0);
        for h_idx in self.mesh.vertex_adjacent_halfedges(v.index, true) {
            let weight = 0.5 * self.dihedral_angle(h_idx) / self.length(self.mesh.halfedges[h_idx].edge.unwrap());
            n.decrement_by(self.vector(h_idx).times(weight));
        }
        n.unit()
    }

    pub fn vertex_normal_mean_curvature(&self, v: &Vertex) -> Vector {
        let mut n = Vector::new(0.0, 0.0, 0.0);
        for h_idx in self.mesh.vertex_adjacent_halfedges(v.index, true) {
            let weight = 0.5 * (self.cotan(h_idx) + self.cotan(self.mesh.halfedges[h_idx].twin.unwrap()));
            n.decrement_by(self.vector(h_idx).times(weight));
        }
        n.unit()
    }

    pub fn vertex_normal_sphere_inscribed(&self, v: &Vertex) -> Vector {
        let mut n = Vector::new(0.0, 0.0, 0.0);
        for c_idx in self.mesh.vertex_adjacent_corners(v.index, true) {
            let c = &self.mesh.corners[c_idx];
            let u = self.vector(self.mesh.halfedges[c.halfedge.unwrap()].prev.unwrap());
            let v_vec = self.vector(self.mesh.halfedges[c.halfedge.unwrap()].next.unwrap()).negated();
            n.increment_by(u.cross(v_vec).over(u.norm2() * v_vec.norm2()));
        }
        n.unit()
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

    pub fn laplace_matrix(&self) -> SparseMatrix {
        let v_count = self.mesh.vertices.len();
        let mut triplet = Triplet::new(v_count, v_count);
        for v in &self.mesh.vertices {
            let i = v.index;
            let mut sum = 0.0;
            for h_idx in self.mesh.vertex_adjacent_halfedges(v.index, true) {
                let j = self.mesh.halfedges[self.mesh.halfedges[h_idx].twin.expect("Twin should exist")].vertex.expect("Vertex should exist");
                let weight = (self.cotan(h_idx) + self.cotan(self.mesh.halfedges[h_idx].twin.expect("Twin should exist"))) / 2.0;
                sum += weight;
                triplet.add_entry(-weight, i, j);
            }
            triplet.add_entry(sum, i, i);
        }
        SparseMatrix::from_triplet(triplet)
    }

    pub fn mass_matrix(&self) -> SparseMatrix {
        let v_count = self.mesh.vertices.len();
        let mut triplet = Triplet::new(v_count, v_count);
        for v in &self.mesh.vertices {
            triplet.add_entry(self.barycentric_dual_area(v), v.index, v.index);
        }
        SparseMatrix::from_triplet(triplet)
    }

    pub fn complex_laplace_matrix(&self) -> ComplexSparseMatrix {
        let v_count = self.mesh.vertices.len();
        let mut triplet = ComplexTriplet::new(v_count, v_count);
        for v in &self.mesh.vertices {
            let i = v.index;
            let mut sum = 0.0;
            for h_idx in self.mesh.vertex_adjacent_halfedges(v.index, true) {
                let j = self.mesh.halfedges[self.mesh.halfedges[h_idx].twin.expect("Twin should exist")].vertex.expect("Vertex should exist");
                let weight = (self.cotan(h_idx) + self.cotan(self.mesh.halfedges[h_idx].twin.expect("Twin should exist"))) / 2.0;
                sum += weight;
                triplet.add_entry(Complex::new(-weight, 0.0), i, j);
            }
            triplet.add_entry(Complex::new(sum, 0.0), i, i);
        }
        ComplexSparseMatrix::from_triplet(triplet)
    }
}
