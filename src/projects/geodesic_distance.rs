use crate::core::geometry::Geometry;
use crate::linear_algebra::{DenseMatrix, SparseMatrix, Vector, Cholesky};
use crate::linear_algebra::traits::{LinearSolver, DenseMatrixOps, Vector3Ops};

pub struct HeatMethod<'a> {
    pub geometry: &'a Geometry<'a>,
    pub a: SparseMatrix<f64>,
    pub f: SparseMatrix<f64>,
}

impl<'a> HeatMethod<'a> {
    pub fn new(geometry: &'a Geometry<'a>) -> Self {
        let t = geometry.mean_edge_length().powi(2);
        let m = geometry.mass_matrix();
        let a = geometry.laplace_matrix();
        let f = &m + &(&a * t);
        Self { geometry, a, f }
    }

    pub fn compute_vector_field(&self, u: &DenseMatrix) -> Vec<Vector> {
        let mut x = vec![faer::Mat::zeros(3, 1); self.geometry.mesh.faces.len()];
        for f in &self.geometry.mesh.faces {
            let normal = self.geometry.face_normal(f).unwrap_or(faer::Mat::zeros(3, 1));
            let area = self.geometry.area(f);
            let mut grad_u = faer::Mat::zeros(3, 1);

            for h_idx in self.geometry.mesh.face_adjacent_halfedges(f.index, true) {
                let prev_idx = self.geometry.mesh.halfedges[h_idx].prev.expect("Halfedge should have a prev");
                let i = self.geometry.mesh.halfedges[prev_idx].vertex.expect("Prev should have a vertex");
                let ui = u[(i, 0)];
                let ei = self.geometry.vector(h_idx);

                let cp = normal.cross(&ei);
                grad_u += &cp * ui;
            }

            grad_u *= 1.0 / (2.0 * area);
            let norm = grad_u.norm();
            if norm > 0.0 {
                grad_u *= 1.0 / norm;
            }
            x[f.index] = -&grad_u;
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
                    let xj = &x[f_idx];
                    let e1 = self.geometry.vector(h_idx);
                    
                    let prev_idx = self.geometry.mesh.halfedges[h_idx].prev.expect("Halfedge should have a prev");
                    let twin_prev_idx = self.geometry.mesh.halfedges[prev_idx].twin.expect("Prev should have a twin");
                    let e2 = self.geometry.vector(twin_prev_idx);
                    
                    let cot_theta1 = self.geometry.cotan(h_idx);
                    let cot_theta2 = self.geometry.cotan(prev_idx);

                    let dot1 = e1.dot(xj);
                    let dot2 = e2.dot(xj);
                    sum += cot_theta1 * dot1 + cot_theta2 * dot2;
                }
            }

            div[(v.index, 0)] = 0.5 * sum;
        }
        div
    }

    fn subtract_minimum_distance(&self, phi: &mut DenseMatrix) {
        let mut min = f64::INFINITY;
        for i in 0..phi.nrows() {
            let val = phi[(i, 0)];
            if val < min { min = val; }
        }

        for i in 0..phi.nrows() {
            phi[(i, 0)] -= min;
        }
    }

    pub fn compute(&self, delta: &DenseMatrix) -> DenseMatrix {
        let llt = Cholesky::<f64>::new(&self.f);
        let u = llt.solve(delta);
        let x = self.compute_vector_field(&u);
        let div = self.compute_divergence(&x);
        
        let llt_a = Cholesky::<f64>::new(&self.a);
        let mut phi = llt_a.solve(&-div);
        self.subtract_minimum_distance(&mut phi);
        phi
    }
}
