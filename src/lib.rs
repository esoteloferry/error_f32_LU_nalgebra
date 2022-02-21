extern crate nalgebra as na;

use std::ops::Mul;

use na::{DMatrix, DVector};
//Test consist in getting eigenvalues and eigenvectors of a 2 dof mass spring system
#[allow(dead_code)]
fn simple_test_lu_f32() {
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

    for (_i, lambda) in eigen.iter().enumerate() {
        let eigen_identity_matrix: DMatrix<f32> = DMatrix::identity(dim, dim).mul(*lambda);
        let a_mat_for_eigen_vec: DMatrix<f32> = a_matrix.clone() - &eigen_identity_matrix;
        let mut b_vec: Vec<f32> = vec![0.0; dim];
        b_vec[1] = 1.0;
        let b: DMatrix<f32> = DMatrix::from_vec(dim, 1, b_vec);
        println!(
            "Matrix used by LU f32 operation: {} which det is : {}",
            a_mat_for_eigen_vec,
            a_mat_for_eigen_vec.determinant()
        );
        let decomp = a_mat_for_eigen_vec.lu();
        match decomp.solve(&b) {
            Some(i) => i,
            None => panic!("Linear resolution failed"),
        };
    }
}
#[allow(dead_code)]
fn simple_test_lu_f64() -> DVector<f64> {
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

    for (_i, lambda) in eigen.iter().enumerate() {
        let eigen_identity_matrix: DMatrix<f64> = DMatrix::identity(dim, dim).mul(*lambda);
        let a_mat_for_eigen_vec: DMatrix<f64> = a_matrix.clone() - &eigen_identity_matrix;
        let mut b_vec: Vec<f64> = vec![0.0; dim];
        b_vec[1] = 1.0;
        let b: DMatrix<f64> = DMatrix::from_vec(dim, 1, b_vec);
        println!(
            "Matrix used by LU f64 operation: {} which det is : {}",
            a_mat_for_eigen_vec,
            a_mat_for_eigen_vec.determinant()
        );
        let decomp = a_mat_for_eigen_vec.lu();
        match decomp.solve(&b) {
            Some(i) => i,
            None => panic!("Linear resolution failed"),
        };
    }

    return eigen;
}

#[allow(dead_code)]
fn simple_test_symm_eigen() -> DVector<f64> {
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

    println!("Matrix used to calculate symmetric_eigen : {}", a_matrix);
    let eigen = a_matrix.symmetric_eigen();
    eigen.eigenvalues
}

#[cfg(test)]
mod tests {
    use crate::*;

    #[test]
    fn lu_f32_test() {
        simple_test_lu_f32();
    }
    #[test]
    fn lu_f64_test() {
        simple_test_lu_f64();
    }
    #[test]
    fn lu_f64_vs_symm_eigen() {
        let lu_f64_eig_values = simple_test_lu_f64();
        let symm_eig_f64_eig_values = simple_test_symm_eigen();

        println!("Eigenvalues from lu f64: {}", lu_f64_eig_values);
        println!(
            "Eigenvalues from symm_eigen f64: {}",
            symm_eig_f64_eig_values
        );
        assert!(lu_f64_eig_values.relative_eq(&symm_eig_f64_eig_values, 0.001, 0.001));
    }
}
