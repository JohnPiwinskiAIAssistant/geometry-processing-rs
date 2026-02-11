use crate::core::geometry::Geometry;
use crate::core::mesh::{MeshBackend, Face};
use crate::linear_algebra::{DenseMatrix, SparseMatrix, Vector, Cholesky};
use crate::linear_algebra::traits::{LinearSolver, DenseMatrixOps, Vector3Ops};

pub struct HeatMethod<'a, B: MeshBackend> {
    pub geometry: &'a Geometry<'a, B>,
    pub a: SparseMatrix<f64>,
    pub f: SparseMatrix<f64>,
}

impl<'a, B: MeshBackend> HeatMethod<'a, B> {
    pub fn new(geometry: &'a Geometry<'a, B>) -> Self {
        let t = geometry.mean_edge_length().powi(2);
        let m = geometry.mass_matrix();
        let a = geometry.laplace_matrix();
        let f = &m + &(&a * t);
        Self { geometry, a, f }
    }

    pub fn compute_vector_field(&self, u: &DenseMatrix) -> Vec<Vector> {
        let num_faces = self.geometry.mesh.num_faces();
        let mut x = vec![faer::Mat::zeros(3, 1); num_faces];
        for f_idx in 0..num_faces {
            let f = Face::new(f_idx);
            let normal = self.geometry.face_normal(&f).unwrap_or(faer::Mat::zeros(3, 1));
            let area = self.geometry.area(&f);
            let mut grad_u = faer::Mat::zeros(3, 1);

            for h_idx in self.geometry.mesh.face_adjacent_halfedges(f_idx, true) {
                let prev_idx = self.geometry.mesh.backend.halfedge_prev(h_idx).expect("Halfedge should have a prev");
                let i = self.geometry.mesh.backend.halfedge_vertex(prev_idx).expect("Prev should have a vertex");
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
        let v_count = self.geometry.mesh.num_vertices();
        let mut div = DenseMatrix::zeros(v_count, 1);

        for i in 0..v_count {
            let mut sum = 0.0;

            for h_idx in self.geometry.mesh.vertex_adjacent_halfedges(i, true) {
                if !self.geometry.mesh.backend.halfedge_on_boundary(h_idx) {
                    let f_idx = self.geometry.mesh.backend.halfedge_face(h_idx).expect("Halfedge should have a face");
                    let xj = &x[f_idx];
                    let e1 = self.geometry.vector(h_idx);
                    
                    let prev_idx = self.geometry.mesh.backend.halfedge_prev(h_idx).expect("Halfedge should have a prev");
                    let twin_prev_idx = self.geometry.mesh.backend.halfedge_twin(prev_idx).expect("Prev should have a twin");
                    let e2 = self.geometry.vector(twin_prev_idx);
                    
                    let cot_theta1 = self.geometry.cotan(h_idx);
                    let cot_theta2 = self.geometry.cotan(prev_idx);

                    let dot1 = e1.dot(xj);
                    let dot2 = e2.dot(xj);
                    sum += cot_theta1 * dot1 + cot_theta2 * dot2;
                }
            }

            div[(i, 0)] = 0.5 * sum;
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
