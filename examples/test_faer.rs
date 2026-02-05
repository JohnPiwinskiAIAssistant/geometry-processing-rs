use faer::sparse::SparseColMat;
use num_complex::Complex64;

fn main() {
    let triplets = vec![
        (0, 0, Complex64::new(1.0, 0.0)),
        (0, 0, Complex64::new(2.0, 0.0)),
    ];
    let mat = SparseColMat::<usize, Complex64>::try_new_from_triplets(
        1, 1, &triplets
    ).unwrap();
    println!("Value: {:?}", mat.values());
}
