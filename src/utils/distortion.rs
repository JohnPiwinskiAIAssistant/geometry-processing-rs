use crate::linear_algebra::Vector;
use crate::core::geometry::Geometry;
use std::collections::HashMap;

pub struct Distortion;

impl Distortion {
    pub fn compute_quasi_conformal_error_per_face(p: &[Vector], q: &[Vector]) -> f64 {
        if p.len() < 3 || q.len() < 3 { return 1.0; }
        
        let u1 = p[1].minus(p[0]);
        let u2 = p[2].minus(p[0]);

        let v1 = q[1].minus(q[0]);
        let v2 = q[2].minus(q[0]);

        let mut e1 = u1;
        e1.normalize();
        let mut e2 = u2.minus(e1.times(u2.dot(e1)));
        e2.normalize();

        let mut f1 = v1;
        f1.normalize();
        let mut f2 = v2.minus(f1.times(v2.dot(f1)));
        f2.normalize();

        let p_proj = [
            Vector::new(0.0, 0.0, 0.0),
            Vector::new(u1.dot(e1), u1.dot(e2), 0.0),
            Vector::new(u2.dot(e1), u2.dot(e2), 0.0),
        ];

        let q_proj = [
            Vector::new(0.0, 0.0, 0.0),
            Vector::new(v1.dot(f1), v1.dot(f2), 0.0),
            Vector::new(v2.dot(f1), v2.dot(f2), 0.0),
        ];

        let a_area = 2.0 * u1.cross(u2).norm();
        if a_area.abs() < 1e-12 { return 1.0; }

        let mut ss = Vector::new(0.0, 0.0, 0.0);
        for i in 0..3 {
            let val = q_proj[i].times(p_proj[(i + 1) % 3].y - p_proj[(i + 2) % 3].y);
            ss.increment_by(val);
        }
        ss.divide_by(a_area);

        let mut st = Vector::new(0.0, 0.0, 0.0);
        for i in 0..3 {
            let val = q_proj[i].times(p_proj[(i + 2) % 3].x - p_proj[(i + 1) % 3].x);
            st.increment_by(val);
        }
        st.divide_by(a_area);

        let a = ss.dot(ss);
        let b = ss.dot(st);
        let c = st.dot(st);
        let det = ((a - c).powi(2) + 4.0 * b * b).sqrt();
        let mut gamma_large = (0.5 * (a + c + det)).sqrt();
        let mut gamma_small = (0.5 * (a + c - det)).sqrt();

        if gamma_large < gamma_small {
            std::mem::swap(&mut gamma_large, &mut gamma_small);
        }

        if gamma_small.abs() < 1e-12 { return 1.0; }
        gamma_large / gamma_small
    }

    pub fn compute_quasi_conformal_error(
        colors: &mut Vec<f64>,
        parameterization: &HashMap<usize, Vector>,
        geometry: &Geometry
    ) -> f64 {
        let mut total_area = 0.0;
        let mut total_qc_error = 0.0;
        let mut qc_errors = HashMap::new();

        for f in &geometry.mesh.faces {
            let mut p = Vec::new();
            let mut q = Vec::new();
            for v_idx in geometry.mesh.face_adjacent_vertices(f.index, true) {
                p.push(geometry.positions[v_idx]);
                q.push(*parameterization.get(&v_idx).unwrap_or(&Vector::new(0.0, 0.0, 0.0)));
            }

            let qc_error = Self::compute_quasi_conformal_error_per_face(&p, &q);
            let area = geometry.area(f);
            total_area += area;
            total_qc_error += qc_error * area;
            qc_errors.insert(f.index, qc_error.max(1.0).min(1.5));
        }

        for v in &geometry.mesh.vertices {
            let i = v.index;
            let mut qc_error = 0.0;
            let mut count = 0;
            for f_idx in geometry.mesh.vertex_adjacent_faces(v.index, true) {
                qc_error += qc_errors.get(&f_idx).cloned().unwrap_or(1.0);
                count += 1;
            }
            if count > 0 {
                qc_error /= count as f64;
            }

            let color = hsv((2.0 - 4.0 * (qc_error - 1.0)) / 3.0, 0.7, 0.65);
            colors[3 * i + 0] = color.x;
            colors[3 * i + 1] = color.y;
            colors[3 * i + 2] = color.z;
        }

        if total_area.abs() < 1e-12 { return 0.0; }
        total_qc_error / total_area
    }

    pub fn compute_area_scaling_per_face(p: &[Vector], q: &[Vector]) -> f64 {
        if p.len() < 3 || q.len() < 3 { return 0.0; }
        let u1 = p[1].minus(p[0]);
        let u2 = p[2].minus(p[0]);
        let big_area = u1.cross(u2).norm();

        let v1 = q[1].minus(q[0]);
        let v2 = q[2].minus(q[0]);
        let small_area = v1.cross(v2).norm();

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
    Vector::new(r, g, b)
}
