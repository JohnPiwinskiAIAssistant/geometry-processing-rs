use crate::core::geometry::Geometry;
use crate::linear_algebra::{DenseMatrix, SparseMatrix, Vector};

pub struct HeatMethod<'a> {
    pub geometry: &'a Geometry<'a>,
    pub a: SparseMatrix,
    pub f: SparseMatrix,
}

impl<'a> HeatMethod<'a> {
    pub fn new(geometry: &'a Geometry<'a>) -> Self {
        let t = geometry.mean_edge_length().powi(2);
        let m = geometry.mass_matrix();
        let a = geometry.laplace_matrix();
        let f = m.plus(&a.times_real(t));
        Self { geometry, a, f }
    }

    pub fn compute_vector_field(&self, u: &DenseMatrix) -> Vec<Vector> {
        let mut x = vec![Vector::new(0.0, 0.0, 0.0); self.geometry.mesh.faces.len()];
        for f in &self.geometry.mesh.faces {
            let normal = self.geometry.face_normal(f).unwrap_or(Vector::new(0.0, 0.0, 0.0));
            let area = self.geometry.area(f);
            let mut grad_u = Vector::new(0.0, 0.0, 0.0);

            for h_idx in self.geometry.mesh.face_adjacent_halfedges(f.index, true) {
                let prev_idx = self.geometry.mesh.halfedges[h_idx].prev.expect("Halfedge should have a prev");
                let i = self.geometry.mesh.halfedges[prev_idx].vertex.expect("Prev should have a vertex");
                let ui = u.get(i, 0);
                let ei = self.geometry.vector(h_idx);

                grad_u.increment_by(normal.cross(ei).times(ui));
            }

            grad_u.divide_by(2.0 * area);
            grad_u.normalize();
            x[f.index] = grad_u.negated();
        }
        x
    }

    pub fn compute_divergence(&self, x: &[Vector]) -> DenseMatrix {
        let v_count = self.geometry.mesh.vertices.len();
        let mut div = DenseMatrix::zeros(v_count, 1);

        for v in &self.geometry.mesh.vertices {
            let mut sum = 0.0;

            for h_idx in self.geometry.mesh.vertex_adjacent_halfedges(v.index, true) {
                if !self.geometry.mesh.halfedges[h_idx].on_boundary {
                    let f_idx = self.geometry.mesh.halfedges[h_idx].face.expect("Halfedge should have a face");
                    let xj = x[f_idx];
                    let e1 = self.geometry.vector(h_idx);
                    
                    let prev_idx = self.geometry.mesh.halfedges[h_idx].prev.expect("Halfedge should have a prev");
                    let twin_prev_idx = self.geometry.mesh.halfedges[prev_idx].twin.expect("Prev should have a twin");
                    let e2 = self.geometry.vector(twin_prev_idx);
                    
                    let cot_theta1 = self.geometry.cotan(h_idx);
                    let cot_theta2 = self.geometry.cotan(prev_idx);

                    sum += cot_theta1 * e1.dot(xj) + cot_theta2 * e2.dot(xj);
                }
            }

            div.set(0.5 * sum, v.index, 0);
        }
        div
    }

    fn subtract_minimum_distance(&self, phi: &mut DenseMatrix) {
        let mut min = f64::INFINITY;
        for i in 0..phi.n_rows() {
            let val = phi.get(i, 0);
            if val < min { min = val; }
        }

        for i in 0..phi.n_rows() {
            phi.set(phi.get(i, 0) - min, i, 0);
        }
    }

    pub fn compute(&self, delta: &DenseMatrix) -> DenseMatrix {
        let u = self.f.chol().solve_positive_definite(delta);
        let x = self.compute_vector_field(&u);
        let div = self.compute_divergence(&x);
        
        let mut phi = self.a.chol().solve_positive_definite(&div.negated());
        self.subtract_minimum_distance(&mut phi);
        phi
    }
}
