use crate::core::geometry::Geometry;
use crate::linear_algebra::{DenseMatrix, SparseMatrix, Cholesky, LU};
use crate::linear_algebra::traits::{SparseOps, LinearSolver, DenseMatrixOps, Vector3Ops};
use crate::projects::vector_field_decomposition::HodgeDecomposition;
use crate::projects::harmonic_bases::HarmonicBases;

pub struct TrivialConnections<'a> {
    pub geometry: &'a Geometry<'a>,
    pub bases: Vec<DenseMatrix>,
    pub p_mat: SparseMatrix,
    pub a_mat: SparseMatrix,
    pub hodge1: SparseMatrix,
    pub d0: SparseMatrix,
}

impl<'a> TrivialConnections<'a> {
    pub fn new(geometry: &'a Geometry<'a>) -> Self {
        let hodge_decomp = HodgeDecomposition::new(geometry);
        let harmonic_bases = HarmonicBases::new(geometry);
        let bases = harmonic_bases.compute(&hodge_decomp);
        
        // Note: Period matrix requires homology generators which should be built by TreeCotree
        // In the JS code, TreeCotree::buildGenerators() is called on the mesh.
        
        let mut tc = Self {
            geometry,
            bases,
            p_mat: crate::linear_algebra::sparse_matrix::identity(0, 0), // stub
            a_mat: hodge_decomp.a,
            hodge1: hodge_decomp.hodge1,
            d0: hodge_decomp.d0,
        };
        tc.p_mat = tc.build_period_matrix();
        tc
    }

    fn build_period_matrix(&self) -> SparseMatrix {
        let n = self.bases.len();
        let mut t = crate::linear_algebra::Triplet::new(n, n);

        for i in 0..n {
            let generator = &self.geometry.mesh.generators[i];
            for j in 0..n {
                let basis = &self.bases[j];
                let mut sum = 0.0;

                for &h_idx in generator {
                    let e_idx = self.geometry.mesh.halfedges[h_idx].edge.expect("Halfedge should have an edge");
                    let sign = if self.geometry.mesh.edges[e_idx].halfedge == Some(h_idx) { 1.0 } else { -1.0 };
                    sum += sign * basis[(e_idx, 0)];
                }
                t.add_entry(sum, i, j);
            }
        }

        SparseOps::from_triplets(n, n, &t.data)
    }

    pub fn satisfy_gauss_bonnet(&self, singularity: &[f64]) -> bool {
        let sum: f64 = singularity.iter().sum();
        (self.geometry.mesh.euler_characteristic() as f64 - sum).abs() < 1e-8
    }

    pub fn compute_co_exact_component(&self, singularity: &[f64]) -> DenseMatrix {
        let v_count = self.geometry.mesh.vertices.len();
        let mut rhs = DenseMatrix::zeros(v_count, 1);
        for i in 0..v_count {
            let u = -self.geometry.angle_defect(&self.geometry.mesh.vertices[i]) + 2.0 * std::f64::consts::PI * singularity[i];
            rhs[(i, 0)] = u;
        }

        let llt = Cholesky::new(&self.a_mat);
        let beta_tilde = llt.solve(&rhs);
        
        &self.hodge1 * &(&self.d0 * &beta_tilde)
    }

    pub fn transport_no_rotation(&self, h_idx: usize, alpha_i: f64) -> f64 {
        let u = self.geometry.vector(h_idx);
        let f_idx = self.geometry.mesh.halfedges[h_idx].face.expect("Should have face");
        let [e1, e2] = self.geometry.orthonormal_bases(&self.geometry.mesh.faces[f_idx]);
        let theta_ij = e2.dot(&u).atan2(e1.dot(&u));

        let twin_idx = self.geometry.mesh.halfedges[h_idx].twin.expect("Should have twin");
        let g_idx = self.geometry.mesh.halfedges[twin_idx].face.expect("Should have face");
        let [f1, f2] = self.geometry.orthonormal_bases(&self.geometry.mesh.faces[g_idx]);
        let theta_ji = f2.dot(&u).atan2(f1.dot(&u));

        alpha_i - theta_ij + theta_ji
    }

    pub fn compute_harmonic_component(&self, delta_beta: &DenseMatrix) -> DenseMatrix {
        let n = self.bases.len();
        let e_count = self.geometry.mesh.edges.len();
        let mut gamma = DenseMatrix::zeros(e_count, 1);
        if n > 0 {
            let mut rhs = DenseMatrix::zeros(n, 1);
            for i in 0..n {
                let generator = &self.geometry.mesh.generators[i];
                let mut sum = 0.0;

                for &h_idx in generator {
                    let e_idx = self.geometry.mesh.halfedges[h_idx].edge.expect("Edge should have an edge");
                    let sign = if self.geometry.mesh.edges[e_idx].halfedge == Some(h_idx) { 1.0 } else { -1.0 };

                    sum += self.transport_no_rotation(h_idx, 0.0);
                    sum -= sign * delta_beta[(e_idx, 0)];
                }

                while sum < -std::f64::consts::PI { sum += 2.0 * std::f64::consts::PI; }
                while sum >= std::f64::consts::PI { sum -= 2.0 * std::f64::consts::PI; }

                rhs[(i, 0)] = sum;
            }

            let lu = LU::new(&self.p_mat);
            let z = lu.solve(&rhs);

            for i in 0..n {
                let basis = &self.bases[i];
                let zi = z[(i, 0)];
                gamma += basis * zi;
            }
        }

        gamma
    }

    pub fn compute_connections(&self, singularity: &[f64]) -> DenseMatrix {
        let delta_beta = self.compute_co_exact_component(singularity);
        let gamma = self.compute_harmonic_component(&delta_beta);
        delta_beta + gamma
    }
}
