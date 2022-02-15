extern crate nalgebra as na;

use std::ops::Mul;

use na::DMatrix;
//Test consist in getting eigenvalues and eigenvectors
fn simple_test_f32() {
    let dim = 2;
    let mass_matrix: DMatrix<f32> = DMatrix::from_vec(dim, dim, vec![1.0, 0.0, 0.0, 2.0]); //kg
    let stiff_matrix: DMatrix<f32> = DMatrix::from_vec(dim, dim, vec![2.0, -1.0, -1.0, 2.0]);
    //N/m
    //

    let inv_m_matrix = match mass_matrix.clone().try_inverse() {
        Some(i) => i,
        None => panic!("Mass matrix could not be inversed"),
    };

    let a_matrix = inv_m_matrix * stiff_matrix;

    let eigen = match a_matrix.eigenvalues() {
        Some(i) => i,
        None => panic!("Error getting eigenvalues"),
    };

    //get first eigenvector

    for (i, lambda) in eigen.iter().enumerate() {
        let eigen_identity_matrix: DMatrix<f32> = DMatrix::identity(dim, dim).mul(*lambda);
        let a_mat_for_eigen_vec: DMatrix<f32> = a_matrix.clone() - &eigen_identity_matrix;
        let mut b_vec: Vec<f32> = vec![0.0; dim];
        b_vec[1] = 1.0;
        let b: DMatrix<f32> = DMatrix::from_vec(dim, 1, b_vec);
        println!("Matrix used by LU operation: {}", a_mat_for_eigen_vec);
        let decomp = a_mat_for_eigen_vec.lu();
        let mut x = match decomp.solve(&b) {
            Some(i) => i,
            None => panic!("Linear resolution failed"),
        };
    }
}

fn simple_test_f64() {
    let dim = 2;
    let mass_matrix: DMatrix<f64> = DMatrix::from_vec(dim, dim, vec![1.0, 0.0, 0.0, 2.0]); //kg
    let stiff_matrix: DMatrix<f64> = DMatrix::from_vec(dim, dim, vec![2.0, -1.0, -1.0, 2.0]);
    //N/m
    //

    let inv_m_matrix = match mass_matrix.clone().try_inverse() {
        Some(i) => i,
        None => panic!("Mass matrix could not be inversed"),
    };

    let a_matrix = inv_m_matrix * stiff_matrix;

    let eigen = match a_matrix.eigenvalues() {
        Some(i) => i,
        None => panic!("Error getting eigenvalues"),
    };

    //get first eigenvector

    for (i, lambda) in eigen.iter().enumerate() {
        let eigen_identity_matrix: DMatrix<f64> = DMatrix::identity(dim, dim).mul(*lambda);
        let a_mat_for_eigen_vec: DMatrix<f64> = a_matrix.clone() - &eigen_identity_matrix;
        let mut b_vec: Vec<f64> = vec![0.0; dim];
        b_vec[1] = 1.0;
        let b: DMatrix<f64> = DMatrix::from_vec(dim, 1, b_vec);
        println!("Matrix used by LU operation: {}", a_mat_for_eigen_vec);
        let decomp = a_mat_for_eigen_vec.lu();
        let mut x = match decomp.solve(&b) {
            Some(i) => i,
            None => panic!("Linear resolution failed"),
        };
    }
}

#[cfg(test)]
mod tests {
    use crate::*;

    #[test]
    fn f32_test() {
        simple_test_f32();
    }
    #[test]
    fn f64_test() {
        simple_test_f64();
    }
}
