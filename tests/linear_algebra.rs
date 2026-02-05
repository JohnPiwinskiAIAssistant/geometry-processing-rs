use geometry_processing_rs::linear_algebra::*;

#[test]
fn test_vector() {
    let u = Vector::new(3.0, 4.0, 0.0);
    assert_eq!(u.x, 3.0);
    assert_eq!(u.y, 4.0);
    assert_eq!(u.z, 0.0);

    let mut v = Vector::new(0.0, 0.0, 0.0);
    v.x = 3.0;
    v.y = 0.0;
    v.z = 3.0;
    assert_eq!(v.x, 3.0);
    assert_eq!(v.y, 0.0);
    assert_eq!(v.z, 3.0);

    assert!((u.norm() - 5.0).abs() < 1e-8);
    assert!((u.norm2() - 25.0).abs() < 1e-8);

    let mut w = u;
    w.normalize();
    assert!((1.0 - w.norm()).abs() < 1e-8);

    let w_unit = u.unit();
    assert!((w_unit.x - 3.0/5.0).abs() < 1e-8);
    assert!((w_unit.y - 4.0/5.0).abs() < 1e-8);
    assert_eq!(w_unit.z, 0.0);

    let inf_v = Vector::new(f64::INFINITY, 3.0, 0.0);
    let nan_v = Vector::new(0.0, 3.0, f64::NAN);
    assert!(!inf_v.is_valid());
    assert!(!nan_v.is_valid());
    assert!(u.is_valid());

    let mut add_v = u;
    add_v.increment_by(v);
    assert_eq!(add_v, Vector::new(6.0, 4.0, 3.0));

    let mut sub_v = u;
    sub_v.decrement_by(v);
    assert_eq!(sub_v, Vector::new(0.0, 4.0, -3.0));

    let mut scale_v = u;
    scale_v.scale_by(2.0);
    assert_eq!(scale_v, Vector::new(6.0, 8.0, 0.0));

    let mut div_v = u;
    div_v.divide_by(2.0);
    assert_eq!(div_v, Vector::new(1.5, 2.0, 0.0));

    assert_eq!(u.plus(v), Vector::new(6.0, 4.0, 3.0));
    assert_eq!(u.minus(v), Vector::new(0.0, 4.0, -3.0));
    assert_eq!(u.times(2.0), Vector::new(6.0, 8.0, 0.0));
    assert_eq!(u.over(2.0), Vector::new(1.5, 2.0, 0.0));
    assert_eq!(u.negated(), Vector::new(-3.0, -4.0, 0.0));

    assert_eq!(u.dot(v), 9.0);
    assert_eq!(u.cross(v), Vector::new(12.0, -9.0, -12.0));
}

#[test]
fn test_dense_matrix() {
    let z = DenseMatrix::zeros(3, 3);
    assert_eq!(z.get(1, 1), 0.0);

    let id = DenseMatrix::identity(3, 3);
    assert_eq!(id.get(0, 0), 1.0);
    assert_eq!(id.get(0, 1), 0.0);

    let ones = DenseMatrix::ones(3, 3);
    assert_eq!(ones.get(2, 2), 1.0);

    let const_mat = DenseMatrix::constant(5.0, 3, 3);
    assert_eq!(const_mat.get(1, 1), 5.0);

    let mut a = DenseMatrix::zeros(2, 2);
    a.set(1.0, 0, 0);
    a.set(2.0, 1, 0);
    a.set(3.0, 0, 1);
    a.set(4.0, 1, 1);
    let at = a.transpose();
    assert_eq!(at.get(0, 1), 2.0);
    assert_eq!(at.get(1, 0), 3.0);

    assert_eq!(a.n_rows(), 2);
    assert_eq!(a.n_cols(), 2);

    let mut n_a = DenseMatrix::zeros(2, 2);
    n_a.set(3.0, 0, 0);
    n_a.set(4.0, 1, 0);
    assert!((n_a.norm(0) - 4.0).abs() < 1e-8); // max row sum
    assert!((n_a.norm(1) - 7.0).abs() < 1e-8); // max col sum
    assert!((n_a.norm(2) - 5.0).abs() < 1e-8); // frobenius

    let rank_id = DenseMatrix::identity(2, 2);
    assert_eq!(rank_id.rank(), 2);

    assert_eq!(a.sum(), 10.0);

    let sub = a.sub_matrix(0, 2, 1, 2);
    assert_eq!(sub.n_rows(), 2);
    assert_eq!(sub.n_cols(), 1);
    assert_eq!(sub.get(0, 0), 3.0);
    assert_eq!(sub.get(1, 0), 4.0);

    let mut scale_a = DenseMatrix::ones(2, 2);
    scale_a.scale_by(5.0);
    assert_eq!(scale_a.get(0, 0), 5.0);

    let mut inc_a = DenseMatrix::ones(2, 2);
    inc_a.increment_by(&DenseMatrix::ones(2, 2));
    assert_eq!(inc_a.get(0, 0), 2.0);

    let hcat = DenseMatrix::ones(2, 2).hcat(&DenseMatrix::ones(2, 2));
    assert_eq!(hcat.n_cols(), 4);
    
    let vcat = DenseMatrix::ones(2, 2).vcat(&DenseMatrix::ones(2, 2));
    assert_eq!(vcat.n_rows(), 4);
}

#[test]
fn test_sparse_matrix() {
    let mut triplet = Triplet::new(2, 2);
    triplet.add_entry(1.0, 0, 0);
    triplet.add_entry(2.0, 1, 0);
    triplet.add_entry(3.0, 0, 1);
    triplet.add_entry(4.0, 1, 1);
    let s = SparseMatrix::from_triplet(triplet);
    let d = s.to_dense();
    assert_eq!(d.get(1, 0), 2.0);
    assert_eq!(d.get(0, 1), 3.0);

    let sid = SparseMatrix::identity(3, 3);
    assert_eq!(sid.get(1, 1), 1.0);
    assert_eq!(sid.nnz(), 3);

    let diag = SparseMatrix::diag(&DenseMatrix::ones(3, 1));
    assert_eq!(diag.get(2, 2), 1.0);

    let st = s.transpose();
    assert_eq!(st.get(0, 1), 2.0);

    let mut triplet_inv = Triplet::new(2, 2);
    triplet_inv.add_entry(1.0, 0, 0);
    triplet_inv.add_entry(2.0, 1, 1);
    let s_inv = SparseMatrix::from_triplet(triplet_inv);
    let sinv = s_inv.invert_diagonal();
    assert_eq!(sinv.get(1, 1), 0.5);

    let mut triplet_norm = Triplet::new(2, 2);
    triplet_norm.add_entry(3.0, 0, 0);
    triplet_norm.add_entry(4.0, 1, 0);
    let s_norm = SparseMatrix::from_triplet(triplet_norm);
    assert!((s_norm.frobenius_norm() - 5.0).abs() < 1e-8);

    // Solvers
    let mut b = DenseMatrix::ones(2, 1);
    b.scale_by(6.0);
    let mut triplet_solv = Triplet::new(2, 2);
    triplet_solv.add_entry(3.0, 0, 0);
    triplet_solv.add_entry(3.0, 1, 1);
    let s_solv = SparseMatrix::from_triplet(triplet_solv);
    
    let x_chol = s_solv.chol().solve_positive_definite(&b);
    assert!((x_chol.get(0, 0) - 2.0).abs() < 1e-8);

    let x_lu = s_solv.lu().solve_square(&b);
    assert!((x_lu.get(1, 0) - 2.0).abs() < 1e-8);

    let x_qr = s_solv.qr().solve(&b);
    assert!((x_qr.get(0, 0) - 2.0).abs() < 1e-8);
}
