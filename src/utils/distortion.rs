use crate::linear_algebra::Vector;
use crate::core::geometry::Geometry;
use crate::core::mesh::{MeshBackend, Face};
use std::collections::HashMap;

pub struct Distortion;

impl Distortion {
    pub fn compute_quasi_conformal_error_per_face(p: &[Vector], q: &[Vector]) -> f64 {
        if p.len() < 3 || q.len() < 3 { return 1.0; }
        
        let u1 = &p[1] - &p[0];
        let u2 = &p[2] - &p[0];

        let v1 = &q[1] - &q[0];
        let v2 = &q[2] - &q[0];

        let mut e1 = u1.clone();
        e1 *= 1.0 / (e1.transpose() * &e1)[(0, 0)].sqrt();
        let mut e2 = &u2 - &(&e1 * (u2.transpose() * &e1)[(0, 0)]);
        e2 *= 1.0 / (e2.transpose() * &e2)[(0, 0)].sqrt();

        let mut f1 = v1.clone();
        f1 *= 1.0 / (f1.transpose() * &f1)[(0, 0)].sqrt();
        let mut f2 = &v2 - &(&f1 * (v2.transpose() * &f1)[(0, 0)]);
        f2 *= 1.0 / (f2.transpose() * &f2)[(0, 0)].sqrt();

        let p_proj = [
            faer::Mat::zeros(3, 1),
            faer::mat![[(u1.transpose() * &e1)[(0, 0)]], [(u1.transpose() * &e2)[(0, 0)]], [0.0]],
            faer::mat![[(u2.transpose() * &e1)[(0, 0)]], [(u2.transpose() * &e2)[(0, 0)]], [0.0]],
        ];

        let q_proj = [
            faer::Mat::zeros(3, 1),
            faer::mat![[(v1.transpose() * &f1)[(0, 0)]], [(v1.transpose() * &f2)[(0, 0)]], [0.0]],
            faer::mat![[(v2.transpose() * &f1)[(0, 0)]], [(v2.transpose() * &f2)[(0, 0)]], [0.0]],
        ];

        let u_cross_v = faer::mat![
            [u1[(1, 0)] * u2[(2, 0)] - u1[(2, 0)] * u2[(1, 0)]],
            [u1[(2, 0)] * u2[(0, 0)] - u1[(0, 0)] * u2[(2, 0)]],
            [u1[(0, 0)] * u2[(1, 0)] - u1[(1, 0)] * u2[(0, 0)]]
        ];
        let a_area: f64 = 2.0 * (u_cross_v.transpose() * &u_cross_v)[(0, 0)].sqrt();
        if a_area.abs() < 1e-12 { return 1.0; }

        let mut ss = faer::Mat::zeros(3, 1);
        for i in 0..3 {
            let val = &q_proj[i] * (p_proj[(i + 1) % 3][(1, 0)] - p_proj[(i + 2) % 3][(1, 0)]);
            ss += &val;
        }
        ss *= 1.0 / a_area;

        let mut st = faer::Mat::zeros(3, 1);
        for i in 0..3 {
            let val = &q_proj[i] * (p_proj[(i + 2) % 3][(0, 0)] - p_proj[(i + 1) % 3][(0, 0)]);
            st += &val;
        }
        st *= 1.0 / a_area;

        let a = (ss.transpose() * &ss)[(0, 0)];
        let b = (ss.transpose() * &st)[(0, 0)];
        let c = (st.transpose() * &st)[(0, 0)];
        let det = ((a - c).powi(2) + 4.0 * b * b).sqrt();
        let mut gamma_large = (0.5 * (a + c + det)).sqrt();
        let mut gamma_small = (0.5 * (a + c - det)).sqrt();

        if gamma_large < gamma_small {
            std::mem::swap(&mut gamma_large, &mut gamma_small);
        }

        if gamma_small.abs() < 1e-12 { return 1.0; }
        gamma_large / gamma_small
    }

    pub fn compute_quasi_conformal_error<B: MeshBackend>(
        colors: &mut Vec<f64>,
        parameterization: &HashMap<usize, Vector>,
        geometry: &Geometry<B>
    ) -> f64 {
        let mut total_area = 0.0;
        let mut total_qc_error = 0.0;
        let mut qc_errors = HashMap::new();

        let num_faces = geometry.mesh.num_faces();
        for f_idx in 0..num_faces {
            let mut p = Vec::new();
            let mut q = Vec::new();
            for v_idx in geometry.mesh.face_adjacent_vertices(f_idx, true) {
                p.push(geometry.positions[v_idx].clone());
                q.push(parameterization.get(&v_idx).cloned().unwrap_or(faer::Mat::zeros(3, 1)));
            }

            let qc_error = Self::compute_quasi_conformal_error_per_face(&p, &q);
            let area = geometry.area(&Face::new(f_idx));
            total_area += area;
            total_qc_error += qc_error * area;
            qc_errors.insert(f_idx, qc_error.max(1.0).min(1.5));
        }

        let num_vertices = geometry.mesh.num_vertices();
        for i in 0..num_vertices {
            let mut qc_error_sum = 0.0;
            let mut count = 0;
            for f_idx in geometry.mesh.vertex_adjacent_faces(i, true) {
                qc_error_sum += qc_errors.get(&f_idx).cloned().unwrap_or(1.0);
                count += 1;
            }
            if count > 0 {
                qc_error_sum /= count as f64;
            }

            let color = hsv((2.0 - 4.0 * (qc_error_sum - 1.0)) / 3.0, 0.7, 0.65);
            colors[3 * i + 0] = color[(0, 0)];
            colors[3 * i + 1] = color[(1, 0)];
            colors[3 * i + 2] = color[(2, 0)];
        }

        if total_area.abs() < 1e-12 { return 0.0; }
        total_qc_error / total_area
    }

    pub fn compute_area_scaling_per_face(p: &[Vector], q: &[Vector]) -> f64 {
        if p.len() < 3 || q.len() < 3 { return 0.0; }
        let u1 = &p[1] - &p[0];
        let u2 = &p[2] - &p[0];
        let big_area_cp = faer::mat![
            [u1[(1, 0)] * u2[(2, 0)] - u1[(2, 0)] * u2[(1, 0)]],
            [u1[(2, 0)] * u2[(0, 0)] - u1[(0, 0)] * u2[(2, 0)]],
            [u1[(0, 0)] * u2[(1, 0)] - u1[(1, 0)] * u2[(0, 0)]]
        ];
        let big_area = (big_area_cp.transpose() * &big_area_cp)[(0, 0)].sqrt();

        let v1 = &q[1] - &q[0];
        let v2 = &q[2] - &q[0];
        let small_area_cp = faer::mat![
            [v1[(1, 0)] * v2[(2, 0)] - v1[(2, 0)] * v2[(1, 0)]],
            [v1[(2, 0)] * v2[(0, 0)] - v1[(0, 0)] * v2[(2, 0)]],
            [v1[(0, 0)] * v2[(1, 0)] - v1[(1, 0)] * v2[(0, 0)]]
        ];
        let small_area = (small_area_cp.transpose() * &small_area_cp)[(0, 0)].sqrt();

        if big_area.abs() < 1e-12 { return 0.0; }
        (small_area / big_area).ln()
    }
}

pub fn hsv(h: f64, s: f64, v: f64) -> Vector {
    let mut r = 0.0;
    let mut g = 0.0;
    let mut b = 0.0;

    if s == 0.0 {
        r = v; g = v; b = v;
    } else {
        let h_prime = if h == 1.0 { 0.0 } else { h } * 6.0;
        let i = h_prime.floor() as i32;
        let f = h_prime - i as f64;
        let p = v * (1.0 - s);
        let q = v * (1.0 - (s * f));
        let t = v * (1.0 - s * (1.0 - f));

        match i {
            0 => { r = v; g = t; b = p; }
            1 => { r = q; g = v; b = p; }
            2 => { r = p; g = v; b = t; }
            3 => { r = p; g = q; b = v; }
            4 => { r = t; g = p; b = v; }
            5 => { r = v; g = p; b = q; }
            _ => {}
        }
    }
    faer::mat![[r], [g], [b]]
}
