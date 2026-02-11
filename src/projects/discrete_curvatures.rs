use crate::core::geometry::Geometry;
use crate::core::mesh::{MeshBackend, Vertex};
use crate::linear_algebra::{DenseMatrix, Vector};

pub struct DiscreteCurvatures<'a, B: MeshBackend> {
    pub geometry: &'a Geometry<'a, B>,
}

impl<'a, B: MeshBackend> DiscreteCurvatures<'a, B> {
    pub fn new(geometry: &'a Geometry<'a, B>) -> Self {
        Self { geometry }
    }

    pub fn compute_gauss_curvatures(&self) -> DenseMatrix {
        let v_count = self.geometry.mesh.num_vertices();
        let mut res = DenseMatrix::zeros(v_count, 1);
        for i in 0..v_count {
            res[(i, 0)] = self.geometry.scalar_gauss_curvature(&Vertex::new(i));
        }
        res
    }

    pub fn compute_mean_curvatures(&self) -> DenseMatrix {
        let v_count = self.geometry.mesh.num_vertices();
        let mut res = DenseMatrix::zeros(v_count, 1);
        for i in 0..v_count {
            res[(i, 0)] = self.geometry.scalar_mean_curvature(&Vertex::new(i));
        }
        res
    }

    pub fn compute_principal_curvatures(&self) -> Vec<[f64; 2]> {
        let num_v = self.geometry.mesh.num_vertices();
        (0..num_v).map(|i| self.geometry.principal_curvatures(&Vertex::new(i))).collect()
    }

    pub fn compute_normals(&self, method: &str) -> Vec<Vector> {
        let num_v = self.geometry.mesh.num_vertices();
        (0..num_v).map(|i| {
            let v = Vertex::new(i);
            match method {
                "Equally Weighted" => self.geometry.vertex_normal_equally_weighted(&v),
                "Area Weighted" => self.geometry.vertex_normal_area_weighted(&v),
                "Angle Weighted" => self.geometry.vertex_normal_angle_weighted(&v),
                "Gauss Curvature" => self.geometry.vertex_normal_gauss_curvature(&v),
                "Mean Curvature" => self.geometry.vertex_normal_mean_curvature(&v),
                "Sphere Inscribed" => self.geometry.vertex_normal_sphere_inscribed(&v),
                _ => faer::Mat::zeros(3, 1),
            }
        }).collect()
    }
}
