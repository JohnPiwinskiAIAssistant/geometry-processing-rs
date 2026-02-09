use geometry_processing_rs::linear_algebra::*;
use geometry_processing_rs::linear_algebra::traits::{SparseOps, LinearSolver, DenseMatrixOps, Vector3Ops};

#[test]
fn test_vector() {
    let u = faer::mat![[3.0], [4.0], [0.0]];
    assert_eq!(u.read(0, 0), 3.0);
    assert_eq!(u.read(1, 0), 4.0);
    assert_eq!(u.read(2, 0), 0.0);

    let mut v = faer::Mat::zeros(3, 1);
    v.write(0, 0, 3.0);
    v.write(1, 0, 0.0);
    v.write(2, 0, 3.0);
    assert_eq!(v.read(0, 0), 3.0);
    assert_eq!(v.read(1, 0), 0.0);
    assert_eq!(v.read(2, 0), 3.0);

    let u_prod: faer::Mat<f64> = u.transpose() * &u;
    let u_norm_sq: f64 = u_prod.read(0, 0);
    let u_norm = u_norm_sq.sqrt();
    assert!((u_norm - 5.0).abs() < 1e-8);
    assert!((u_norm_sq - 25.0).abs() < 1e-8);

    let mut w = u.clone();
    let w_norm_prod: faer::Mat<f64> = w.transpose() * &w;
    let w_norm_sq_val: f64 = w_norm_prod.read(0, 0);
    let w_norm = w_norm_sq_val.sqrt();
    w *= 1.0 / w_norm;
    let w_unit_prod: faer::Mat<f64> = w.transpose() * &w;
    let w_unit_norm_sq_val: f64 = w_unit_prod.read(0, 0);
    let w_unit_norm = w_unit_norm_sq_val.sqrt();
    assert!((1.0 - w_unit_norm).abs() < 1e-8);

    let u_norm_prod: faer::Mat<f64> = u.transpose() * &u;
    let u_norm_val: f64 = u_norm_prod.read(0, 0).sqrt();
    let w_unit: Vector = &u * (1.0 / u_norm_val);
    assert!((w_unit.read(0, 0) - 3.0/5.0).abs() < 1e-8);
    assert!((w_unit.read(1, 0) - 4.0/5.0).abs() < 1e-8);
    assert_eq!(w_unit.read(2, 0), 0.0);

    // is_valid check
    let inf_v = faer::mat![[f64::INFINITY], [3.0], [0.0]];
    let nan_v = faer::mat![[0.0], [3.0], [f64::NAN]];
    let is_valid = |m: &faer::Mat<f64>| {
        for j in 0..m.ncols() {
            for i in 0..m.nrows() {
                if !m.read(i, j).is_finite() { return false; }
            }
        }
        true
    };
    assert!(!is_valid(&inf_v));
    assert!(!is_valid(&nan_v));
    assert!(is_valid(&u));

    let mut add_v = u.clone();
    add_v += &v;
    assert_eq!(add_v, faer::mat![[6.0], [4.0], [3.0]]);

    let mut sub_v = u.clone();
    sub_v -= &v;
    assert_eq!(sub_v, faer::mat![[0.0], [4.0], [-3.0]]);

    let mut scale_v = u.clone();
    scale_v *= 2.0;
    assert_eq!(scale_v, faer::mat![[6.0], [8.0], [0.0]]);

    let mut div_v = u.clone();
    div_v *= 0.5;
    assert_eq!(div_v, faer::mat![[1.5], [2.0], [0.0]]);

    assert_eq!(&u + &v, faer::mat![[6.0], [4.0], [3.0]]);
    assert_eq!(&u - &v, faer::mat![[0.0], [4.0], [-3.0]]);
    assert_eq!(&u * 2.0, faer::mat![[6.0], [8.0], [0.0]]);
    assert_eq!(&u * 0.5, faer::mat![[1.5], [2.0], [0.0]]);
    assert_eq!(-&u, faer::mat![[-3.0], [-4.0], [0.0]]);

    let dot_prod: faer::Mat<f64> = u.transpose() * &v;
    let dot: f64 = dot_prod.read(0, 0);
    assert_eq!(dot, 9.0);
    
    let cross = faer::mat![
        [u.read(1, 0) * v.read(2, 0) - u.read(2, 0) * v.read(1, 0)],
        [u.read(2, 0) * v.read(0, 0) - u.read(0, 0) * v.read(2, 0)],
        [u.read(0, 0) * v.read(1, 0) - u.read(1, 0) * v.read(0, 0)]
    ];
    assert_eq!(cross, faer::mat![[12.0], [-9.0], [-12.0]]);
}

#[test]
fn test_dense_matrix() {
    let z = DenseMatrix::zeros(3, 3);
    assert_eq!(z[(1, 1)], 0.0);

    let mut id = DenseMatrix::zeros(3, 3);
    for i in 0..3 { id[(i, i)] = 1.0; }
    assert_eq!(id[(0, 0)], 1.0);
    assert_eq!(id[(0, 1)], 0.0);

    let ones = DenseMatrix::from_fn(3, 3, |_, _| 1.0);
    assert_eq!(ones[(2, 2)], 1.0);

    let const_mat = DenseMatrix::from_fn(3, 3, |_, _| 5.0);
    assert_eq!(const_mat[(1, 1)], 5.0);

    let mut a = DenseMatrix::zeros(2, 2);
    a[(0, 0)] = 1.0;
    a[(1, 0)] = 2.0;
    a[(0, 1)] = 3.0;
    a[(1, 1)] = 4.0;
    let at = a.transpose().to_owned();
    assert_eq!(at[(0, 1)], 2.0);
    assert_eq!(at[(1, 0)], 3.0);

    assert_eq!(a.nrows(), 2);
    assert_eq!(a.ncols(), 2);

    let mut n_a = DenseMatrix::zeros(2, 2);
    n_a[(0, 0)] = 3.0;
    n_a[(1, 0)] = 4.0;
    
    let norm_max = (0..n_a.nrows()).map(|i| (0..n_a.ncols()).map(|j| n_a[(i, j)].abs()).sum::<f64>()).fold(0.0, f64::max);
    let norm_l1 = (0..n_a.ncols()).map(|j| (0..n_a.nrows()).map(|i| n_a[(i, j)].abs()).sum::<f64>()).fold(0.0, f64::max);
    assert!((norm_max - 4.0).abs() < 1e-8); 
    assert!((norm_l1 - 7.0).abs() < 1e-8); 
    
    // Manual Frobenius norm for faer::Mat
    let mut frob_sq = 0.0;
    for j in 0..n_a.ncols() {
        for i in 0..n_a.nrows() {
            frob_sq += n_a[(i, j)] * n_a[(i, j)];
        }
    }
    assert!((frob_sq.sqrt() - 5.0).abs() < 1e-8);

    let mut sum = 0.0;
    for j in 0..a.ncols() {
        for i in 0..a.nrows() {
            sum += a[(i, j)];
        }
    }
    assert_eq!(sum, 10.0);

    let sub = a.submatrix(0, 0, 2, 1).to_owned();
    assert_eq!(sub.nrows(), 2);
    assert_eq!(sub.ncols(), 1);
    assert_eq!(sub[(0, 0)], 1.0);
    assert_eq!(sub[(1, 0)], 2.0);

    let mut scale_a = DenseMatrix::from_fn(2, 2, |_, _| 1.0);
    scale_a *= 5.0;
    assert_eq!(scale_a[(0, 0)], 5.0);

    let mut inc_a = DenseMatrix::from_fn(2, 2, |_, _| 1.0);
    inc_a += &DenseMatrix::from_fn(2, 2, |_, _| 1.0);
    assert_eq!(inc_a[(0, 0)], 2.0);

    let hcat = faer::concat![[DenseMatrix::from_fn(2, 2, |_, _| 1.0), DenseMatrix::from_fn(2, 2, |_, _| 1.0)]];
    assert_eq!(hcat.ncols(), 4);
    
    let vcat = faer::concat![[DenseMatrix::from_fn(2, 2, |_, _| 1.0)], [DenseMatrix::from_fn(2, 2, |_, _| 1.0)]];
    assert_eq!(vcat.nrows(), 4);
}

#[test]
fn test_sparse_matrix() {
    use geometry_processing_rs::linear_algebra::{DenseMatrix, SparseMatrix};
    use geometry_processing_rs::linear_algebra::sparse_matrix::{identity, diag, Triplet, Cholesky, LU, QR};
    let mut triplet = Triplet::<f64>::new(2, 2);
    triplet.add_entry(1.0, 0, 0);
    triplet.add_entry(2.0, 1, 0);
    triplet.add_entry(3.0, 0, 1);
    triplet.add_entry(4.0, 1, 1);
    let s = SparseMatrix::<f64>::from_triplets(2, 2, &triplet.data);
    let d = s.to_dense();
    assert_eq!(d[(1, 0)], 2.0);
    assert_eq!(d[(0, 1)], 3.0);

    let sid = identity::<f64>(3, 3);
    assert_eq!(sid.get_val(1, 1), 1.0);
    assert_eq!(sid.compute_nnz(), 3);

    let diag_mat = diag(&DenseMatrix::from_fn(3, 1, |_, _| 1.0));
    assert_eq!(diag_mat.get_val(2, 2), 1.0);

    let st = s.transpose().to_col_major().unwrap();
    assert_eq!(st.get_val(0, 1), 2.0);

    let mut triplet_inv = Triplet::<f64>::new(2, 2);
    triplet_inv.add_entry(1.0, 0, 0);
    triplet_inv.add_entry(2.0, 1, 1);
    let s_inv = SparseMatrix::<f64>::from_triplets(2, 2, &triplet_inv.data);
    let sinv = s_inv.invert_diagonal();
    assert_eq!(sinv.get_val(1, 1), 0.5);

    let mut triplet_norm = Triplet::<f64>::new(2, 2);
    triplet_norm.add_entry(3.0, 0, 0);
    triplet_norm.add_entry(4.0, 1, 0);
    let s_norm = SparseMatrix::<f64>::from_triplets(2, 2, &triplet_norm.data);
    assert!((s_norm.frobenius_norm() - 5.0).abs() < 1e-8);

    // Solvers
    let b = DenseMatrix::from_fn(2, 1, |_, _| 6.0);
    let mut triplet_solv = Triplet::<f64>::new(2, 2);
    triplet_solv.add_entry(3.0, 0, 0);
    triplet_solv.add_entry(3.0, 1, 1);
    let s_solv = SparseMatrix::<f64>::from_triplets(2, 2, &triplet_solv.data);
    
    let x_chol = Cholesky::<f64>::new(&s_solv).solve(&b);
    assert!((x_chol.read(0, 0) - 2.0).abs() < 1e-8);
    
    let x_lu = LU::<f64>::new(&s_solv).solve(&b);
    assert!((x_lu.read(1, 0) - 2.0).abs() < 1e-8);
    
    let x_qr = QR::<f64>::new(&s_solv).solve(&b);
    assert!((x_qr.read(0, 0) - 2.0).abs() < 1e-8);
}
