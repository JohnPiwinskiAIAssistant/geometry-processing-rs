use crate::linear_algebra::Vector;
use crate::utils::distortion::hsv;

pub struct Colormap;

impl Colormap {
    pub fn colormap(x: f64, min: f64, max: f64, values: &[(f64, [f64; 3])]) -> Vector {
        let mut x = x.max(min).min(max);
        if max > min {
            x = (x - min) / (max - min);
        } else {
            x = 0.0;
        }

        let mut i = 1;
        while i < values.len() && values[i].0 < x {
            i += 1;
        }
        if i >= values.len() { i = values.len() - 1; }
        let i = i - 1;

        let c1 = Vector::new(values[i].1[0], values[i].1[1], values[i].1[2]);
        let c2 = Vector::new(values[i+1].1[0], values[i+1].1[1], values[i+1].1[2]);
        
        let denom = (values[i].0 - values[i+1].0).abs();
        let scaling = if denom > 1e-12 { (x - values[i].0) / denom } else { 0.0 };

        c1.plus(c2.minus(c1).times(scaling))
    }

    pub fn seismic_data() -> Vec<(f64, [f64; 3])> {
        vec![
            (0.000, [0.000, 0.000, 0.300]),
            (0.250, [0.004, 0.004, 1.000]),
            (0.500, [1.000, 0.992, 0.992]),
            (0.750, [1.000, 0.004, 0.004]),
            (1.000, [0.500, 0.000, 0.000]),
        ]
    }
}
