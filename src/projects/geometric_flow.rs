use crate::core::geometry::Geometry;
use crate::linear_algebra::{DenseMatrix, SparseMatrix};

pub struct MeanCurvatureFlow<'a, 'b> {
    pub geometry: &'a mut Geometry<'b>,
}

impl<'a, 'b> MeanCurvatureFlow<'a, 'b> {
    pub fn new(geometry: &'a mut Geometry<'b>) -> Self {
        Self { geometry }
    }

    pub fn build_flow_operator(&self, m: &SparseMatrix, h: f64) -> SparseMatrix {
        let a = self.geometry.laplace_matrix();
        // F = M + hA
        m.plus(&a.times_real(h))
    }

    pub fn integrate(&mut self, h: f64) {
        let v_count = self.geometry.mesh.vertices.len();
        let m = self.geometry.mass_matrix();
        let f = self.build_flow_operator(&m, h);

        // construct right hand side
        let mut f0 = DenseMatrix::zeros(v_count, 3);
        for v in &self.geometry.mesh.vertices {
            let i = v.index;
            let p = self.geometry.positions[i];

            f0.set(p.x, i, 0);
            f0.set(p.y, i, 1);
            f0.set(p.z, i, 2);
        }

        let rhs = m.times_dense(&f0);

        // solve linear system (M + hA)fh = Mf0
        let llt = f.chol();
        let fh = llt.solve_positive_definite(&rhs);

        // update positions
        for v in &self.geometry.mesh.vertices {
            let i = v.index;
            self.geometry.positions[i].x = fh.get(i, 0);
            self.geometry.positions[i].y = fh.get(i, 1);
            self.geometry.positions[i].z = fh.get(i, 2);
        }

        // center mesh positions around origin
        self.geometry.normalize(false);
    }
}

pub struct ModifiedMeanCurvatureFlow<'a, 'b> {
    pub geometry: &'a mut Geometry<'b>,
    pub laplace: SparseMatrix,
}

impl<'a, 'b> ModifiedMeanCurvatureFlow<'a, 'b> {
    pub fn new(geometry: &'a mut Geometry<'b>) -> Self {
        let laplace = geometry.laplace_matrix();
        Self { geometry, laplace }
    }

    pub fn build_flow_operator(&self, m: &SparseMatrix, h: f64) -> SparseMatrix {
        // F = M + hA
        m.plus(&self.laplace.times_real(h))
    }

    pub fn integrate(&mut self, h: f64) {
        let v_count = self.geometry.mesh.vertices.len();
        let m = self.geometry.mass_matrix();
        let f = self.build_flow_operator(&m, h);

        let mut f0 = DenseMatrix::zeros(v_count, 3);
        for v in &self.geometry.mesh.vertices {
            let i = v.index;
            let p = self.geometry.positions[i];

            f0.set(p.x, i, 0);
            f0.set(p.y, i, 1);
            f0.set(p.z, i, 2);
        }

        let rhs = m.times_dense(&f0);

        let llt = f.chol();
        let fh = llt.solve_positive_definite(&rhs);

        for v in &self.geometry.mesh.vertices {
            let i = v.index;
            self.geometry.positions[i].x = fh.get(i, 0);
            self.geometry.positions[i].y = fh.get(i, 1);
            self.geometry.positions[i].z = fh.get(i, 2);
        }

        self.geometry.normalize(false);
    }
}
