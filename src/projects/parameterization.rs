use crate::core::geometry::Geometry;
use crate::linear_algebra::{DenseMatrix, SparseMatrix, ComplexSparseMatrix, Vector, ComplexTriplet, Cholesky};
use crate::linear_algebra::traits::{Scalar, SparseOps, LinearSolver, DenseMatrixOps};
use num_complex::Complex64;
use crate::utils::solvers::Solvers;
use std::collections::HashMap;

pub struct SpectralConformalParameterization<'a> {
    pub geometry: &'a Geometry<'a>,
}

impl<'a> SpectralConformalParameterization<'a> {
    pub fn new(geometry: &'a Geometry<'a>) -> Self {
        Self { geometry }
    }

    pub fn build_conformal_energy(&self) -> ComplexSparseMatrix {
        let ed = self.geometry.complex_laplace_matrix();
        let ed = ed.scale(Complex64::new(0.5, 0.0));

        let i_img = Complex64::new(0.0, 1.0);
        let n = ed.nrows();
        let mut t = ComplexTriplet::new(n, n);
        for b in &self.geometry.mesh.boundaries {
            let start = b.halfedge.expect("Boundary should have halfedge");
            let mut curr = start;
            loop {
                let h_idx = curr;
                let v_idx = self.geometry.mesh.halfedges[h_idx].vertex.expect("Should have vertex");
                let twin_idx = self.geometry.mesh.halfedges[h_idx].twin.expect("Should have twin");
                let twin_v_idx = self.geometry.mesh.halfedges[twin_idx].vertex.expect("Should have vertex");

                t.add_entry(i_img * 0.25, v_idx, twin_v_idx);
                t.add_entry(i_img * -0.25, twin_v_idx, v_idx);

                curr = self.geometry.mesh.halfedges[curr].next.expect("Next should exist");
                if curr == start { break; }
            }
        }

        let a = ComplexSparseMatrix::from_triplets(n, n, &t.data);
        &ed - &a
    }

    pub fn flatten(&self) -> Vec<Vector> {
        let ec = self.build_conformal_energy();
        let z = Solvers::solve_inverse_power_method(&ec);
        
        let mut flattening = vec![faer::Mat::zeros(3, 1); self.geometry.mesh.vertices.len()];
        for i in 0..self.geometry.mesh.vertices.len() {
            let zi = z.read(i, 0);
            flattening[i] = faer::mat![[zi.re], [zi.im], [0.0]];
        }

        self.normalize_flattening(&mut flattening);
        flattening
    }

    fn normalize_flattening(&self, flattening: &mut [Vector]) {
        let n = flattening.len();
        if n == 0 { return; }
        
        let mut cm = faer::Mat::zeros(3, 1);
        for p in flattening.iter() {
            cm += p;
        }
        cm *= 1.0 / (n as f64);

        let mut radius: f64 = -1.0;
        for p in flattening.iter_mut() {
            *p -= &cm;
            radius = radius.max((p.transpose() * &*p)[(0, 0)].sqrt());
        }

        if radius > 0.0 {
            for p in flattening.iter_mut() {
                *p *= 1.0 / radius;
            }
        }
    }
}

pub struct BoundaryFirstFlattening<'a> {
    pub geometry: &'a Geometry<'a>,
    pub n_v: usize,
    pub n_i: usize,
    pub n_b: usize,
    pub vertex_index: HashMap<usize, usize>,
    pub b_vertex_index: HashMap<usize, usize>,
    pub k_gaussian: DenseMatrix,
    pub k_geodesic: DenseMatrix,
    pub l_lengths: DenseMatrix,
    pub a_ii: SparseMatrix,
    pub a_ib: SparseMatrix,
    pub a_bb: SparseMatrix,
    pub a_full: SparseMatrix,
}

impl<'a> BoundaryFirstFlattening<'a> {
    pub fn new(geometry: &'a Geometry<'a>) -> Self {
        let mut bff = Self {
            geometry,
            n_v: 0, n_i: 0, n_b: 0,
            vertex_index: HashMap::new(),
            b_vertex_index: HashMap::new(),
            k_gaussian: DenseMatrix::zeros(0, 0),
            k_geodesic: DenseMatrix::zeros(0, 0),
            l_lengths: DenseMatrix::zeros(0, 0),
            a_ii: crate::linear_algebra::sparse_matrix::identity(0, 0),
            a_ib: crate::linear_algebra::sparse_matrix::identity(0, 0),
            a_bb: crate::linear_algebra::sparse_matrix::identity(0, 0),
            a_full: crate::linear_algebra::sparse_matrix::identity(0, 0),
        };

        bff.index_vertices();
        bff.compute_integrated_curvatures();
        bff.compute_boundary_lengths();

        let a = bff.build_special_laplace();
        bff.a_full = a;
        
        bff.a_ii = bff.a_full.sub_matrix(0, bff.n_i, 0, bff.n_i);
        bff.a_ib = bff.a_full.sub_matrix(0, bff.n_i, bff.n_i, bff.n_v);
        bff.a_bb = bff.a_full.sub_matrix(bff.n_i, bff.n_v, bff.n_i, bff.n_v);

        bff
    }

    fn index_vertices(&mut self) {
        self.n_v = self.geometry.mesh.vertices.len();
        for v in &self.geometry.mesh.vertices {
            if !self.geometry.mesh.on_boundary(v.index) {
                self.vertex_index.insert(v.index, self.n_i);
                self.n_i += 1;
            }
        }
        for v in &self.geometry.mesh.vertices {
            if self.geometry.mesh.on_boundary(v.index) {
                self.b_vertex_index.insert(v.index, self.n_b);
                self.vertex_index.insert(v.index, self.n_i + self.n_b);
                self.n_b += 1;
            }
        }
    }

    fn compute_integrated_curvatures(&mut self) {
        self.k_gaussian = DenseMatrix::zeros(self.n_i, 1);
        self.k_geodesic = DenseMatrix::zeros(self.n_b, 1);
        for v in &self.geometry.mesh.vertices {
            let ad = self.geometry.angle_defect(v);
            if self.geometry.mesh.on_boundary(v.index) {
                if let Some(&bi) = self.b_vertex_index.get(&v.index) {
                    self.k_geodesic[(bi, 0)] = ad;
                }
            } else {
                if let Some(&i) = self.vertex_index.get(&v.index) {
                    self.k_gaussian[(i, 0)] = ad;
                }
            }
        }
    }

    fn compute_boundary_lengths(&mut self) {
        self.l_lengths = DenseMatrix::zeros(self.n_b, 1);
        if self.geometry.mesh.boundaries.is_empty() { return; }
        let b = &self.geometry.mesh.boundaries[0];
        let start = b.halfedge.expect("Boundary should have halfedge");
        let mut curr = start;
        loop {
            let h_idx = curr;
            let v_idx = self.geometry.mesh.halfedges[h_idx].vertex.expect("Vertex should exist");
            if let Some(&bi) = self.b_vertex_index.get(&v_idx) {
                let e_idx = self.geometry.mesh.halfedges[h_idx].edge.expect("Edge should exist");
                self.l_lengths[(bi, 0)] = self.geometry.length(e_idx);
            }
            curr = self.geometry.mesh.halfedges[curr].next.expect("Next should exist");
            if curr == start { break; }
        }
    }

    fn build_special_laplace(&self) -> SparseMatrix {
        let mut t = crate::linear_algebra::Triplet::new(self.n_v, self.n_v);
        for v in &self.geometry.mesh.vertices {
            let i = *self.vertex_index.get(&v.index).unwrap();
            let mut sum = 1e-8;
            for h_idx in self.geometry.mesh.vertex_adjacent_halfedges(v.index, true) {
                let twin_idx = self.geometry.mesh.halfedges[h_idx].twin.unwrap();
                let j_orig = self.geometry.mesh.halfedges[twin_idx].vertex.unwrap();
                let j = *self.vertex_index.get(&j_orig).unwrap();
                let weight = (self.geometry.cotan(h_idx) + self.geometry.cotan(twin_idx)) / 2.0;
                sum += weight;
                t.add_entry(-weight, i, j);
            }
            t.add_entry(sum, i, i);
        }
        SparseOps::from_triplets(self.n_v, self.n_v, &t.data)
    }

    pub fn flatten(&self, target: &DenseMatrix, given_scale_factors: bool) -> Vec<Vector> {
        let (u, ktilde) = if given_scale_factors {
            let u = target.clone();
            let h = self.dirichlet_to_neumann(&-target, &u);
            let ktilde = &self.k_geodesic - h;
            (u, ktilde)
        } else {
            let ktilde = target.clone();
            let h = &self.k_geodesic - &ktilde;
            let u = self.neumann_to_dirichlet(&-&self.k_gaussian, &h);
            (u, ktilde)
        };

        self.flatten_with_sf_and_curvatures(&u, &ktilde)
    }

    fn dirichlet_to_neumann(&self, phi: &DenseMatrix, g: &DenseMatrix) -> DenseMatrix {
        let llt = Cholesky::new(&self.a_ii);
        let rhs = phi - &(&self.a_ib * g);
        let a = llt.solve(&rhs);
        let first = self.a_ib.transpose().to_col_major().unwrap() * &a;
        let second = &self.a_bb * g;
        -(first + second)
    }

    fn neumann_to_dirichlet(&self, phi: &DenseMatrix, h: &DenseMatrix) -> DenseMatrix {
        let llt = Cholesky::new(&self.a_full);
        let rhs = faer::concat![[phi], [-h]];
        let a = llt.solve(&rhs);
        a.submatrix(self.n_i, 0, self.n_v - self.n_i, 1).to_owned()
    }

    fn flatten_with_sf_and_curvatures(&self, u: &DenseMatrix, ktilde: &DenseMatrix) -> Vec<Vector> {
        let _u = u;
        let _kt = ktilde;
        vec![faer::Mat::zeros(3, 1); self.n_v]
    }
}
