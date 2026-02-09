use crate::core::geometry::Geometry;
use crate::linear_algebra::{DenseMatrix, Vector};

pub struct DiscreteCurvatures<'a> {
    pub geometry: &'a Geometry<'a>,
}

impl<'a> DiscreteCurvatures<'a> {
    pub fn new(geometry: &'a Geometry<'a>) -> Self {
        Self { geometry }
    }

    pub fn compute_gauss_curvatures(&self) -> DenseMatrix {
        let v_count = self.geometry.mesh.vertices.len();
        let mut res = DenseMatrix::zeros(v_count, 1);
        for v in &self.geometry.mesh.vertices {
            res[(v.index, 0)] = self.geometry.scalar_gauss_curvature(v);
        }
        res
    }

    pub fn compute_mean_curvatures(&self) -> DenseMatrix {
        let v_count = self.geometry.mesh.vertices.len();
        let mut res = DenseMatrix::zeros(v_count, 1);
        for v in &self.geometry.mesh.vertices {
            res[(v.index, 0)] = self.geometry.scalar_mean_curvature(v);
        }
        res
    }

    pub fn compute_principal_curvatures(&self) -> Vec<[f64; 2]> {
        self.geometry.mesh.vertices.iter()
            .map(|v| self.geometry.principal_curvatures(v))
            .collect()
    }

    pub fn compute_normals(&self, method: &str) -> Vec<Vector> {
        self.geometry.mesh.vertices.iter().map(|v| {
            match method {
                "Equally Weighted" => self.geometry.vertex_normal_equally_weighted(v),
                "Area Weighted" => self.geometry.vertex_normal_area_weighted(v),
                "Angle Weighted" => self.geometry.vertex_normal_angle_weighted(v),
                "Gauss Curvature" => self.geometry.vertex_normal_gauss_curvature(v),
                "Mean Curvature" => self.geometry.vertex_normal_mean_curvature(v),
                "Sphere Inscribed" => self.geometry.vertex_normal_sphere_inscribed(v),
                _ => faer::Mat::zeros(3, 1),
            }
        }).collect()
    }
}
