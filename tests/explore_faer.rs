use faer::sparse::SparseColMat;
use num_complex::Complex64;

#[test]
fn test_faer_sparse_capabilities() {
    let tri = vec![(0usize, 0usize, 2.0f64)];
    let mat = SparseColMat::<usize, f64>::try_new_from_triplets(1, 1, &tri).unwrap();
    
    // Test if addition works
    let sum = &mat + &mat; 
    println!("Addition values: {:?}", sum.values());
    
    // Test if multiplication works
    let prod = &mat * &mat;
    println!("Multiplication values: {:?}", prod.values());
    
    // Test if transpose works
    let trans = mat.transpose();
    println!("Transpose NNZ: {}", trans.compute_nnz());
}
