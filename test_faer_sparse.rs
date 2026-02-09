use faer::sparse::SparseColMat;
use num_complex::Complex64;

fn main() {
    let tri = vec![(0, 0, 1.0), (0, 0, 2.0)];
    let mat = SparseColMat::try_new_from_triplets(1, 1, &tri).unwrap();
    println!("NNZ: {}", mat.compute_nnz());
    // let result = &mat + &mat; // Test if addition works
}
