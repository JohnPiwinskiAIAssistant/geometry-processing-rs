use crate::core::geometry::Geometry;
use crate::core::mesh::MeshBackend;
use crate::linear_algebra::{DenseMatrix, SparseMatrix, Cholesky};
use crate::linear_algebra::traits::{SparseOps, LinearSolver};

pub struct MeanCurvatureFlow;

impl MeanCurvatureFlow {
    pub fn new() -> Self {
        Self
    }

    pub fn build_flow_operator<B: MeshBackend>(&self, geometry: &Geometry<B>, m: &SparseMatrix<f64>, h: f64) -> SparseMatrix<f64> {
        let a = geometry.laplace_matrix();
        m + &a.scale(h)
    }

    pub fn integrate<B: MeshBackend>(&self, geometry: &mut Geometry<B>, h: f64) {
        let v_count = geometry.mesh.num_vertices();
        let m = geometry.mass_matrix();
        let f = self.build_flow_operator(geometry, &m, h);

        // construct right hand side
        let mut f0 = DenseMatrix::zeros(v_count, 3);
        for i in 0..v_count {
            let p = &geometry.positions[i];

            f0[(i, 0)] = p[(0, 0)];
            f0[(i, 1)] = p[(1, 0)];
            f0[(i, 2)] = p[(2, 0)];
        }

        let rhs = &m * &f0;

        // solve linear system (M + hA)fh = Mf0
        let llt = Cholesky::<f64>::new(&f);
        let fh = llt.solve(&rhs);

        // update positions
        for i in 0..v_count {
            geometry.positions[i][(0, 0)] = fh[(i, 0)];
            geometry.positions[i][(1, 0)] = fh[(i, 1)];
            geometry.positions[i][(2, 0)] = fh[(i, 2)];
        }

        // center mesh positions around origin
        geometry.normalize(false);
    }
}

pub struct ModifiedMeanCurvatureFlow {
    pub laplace: SparseMatrix<f64>,
}

impl ModifiedMeanCurvatureFlow {
    pub fn new<B: MeshBackend>(geometry: &Geometry<B>) -> Self {
        let laplace = geometry.laplace_matrix();
        Self { laplace }
    }

    pub fn build_flow_operator(&self, m: &SparseMatrix<f64>, h: f64) -> SparseMatrix<f64> {
        m + &self.laplace.scale(h)
    }

    pub fn integrate<B: MeshBackend>(&self, geometry: &mut Geometry<B>, h: f64) {
        let v_count = geometry.mesh.num_vertices();
        let m = geometry.mass_matrix();
        let f = self.build_flow_operator(&m, h);

        let mut f0 = DenseMatrix::zeros(v_count, 3);
        for i in 0..v_count {
            let p = &geometry.positions[i];

            f0[(i, 0)] = p[(0, 0)];
            f0[(i, 1)] = p[(1, 0)];
            f0[(i, 2)] = p[(2, 0)];
        }

        let rhs = &m * &f0;

        let llt = Cholesky::<f64>::new(&f);
        let fh = llt.solve(&rhs);

        for i in 0..v_count {
            geometry.positions[i][(0, 0)] = fh[(i, 0)];
            geometry.positions[i][(1, 0)] = fh[(i, 1)];
            geometry.positions[i][(2, 0)] = fh[(i, 2)];
        }

        geometry.normalize(false);
    }
}
