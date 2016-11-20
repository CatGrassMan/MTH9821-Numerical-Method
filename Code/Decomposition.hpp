#ifndef DECOMPOSITION_HPP
#define DECOMPOSITION_HPP

#include<Eigen/Dense>
#include<tuple>


//Matrix Multiplication
Eigen::MatrixXd matrix_multiply(const Eigen::MatrixXd &A, const Eigen::MatrixXd &B);

//LU Decomposition without pivoting
std::tuple<Eigen::MatrixXd,Eigen::MatrixXd> lu_no_pivoting(Eigen::MatrixXd &A);

//LU Decomposition with row pivoting
std::tuple<Eigen::MatrixXd,Eigen::MatrixXd, Eigen::MatrixXd> lu_row_pivoting(Eigen::MatrixXd &A);

//Cholesky
Eigen::MatrixXd cholesky(Eigen::MatrixXd &A);


#endif






