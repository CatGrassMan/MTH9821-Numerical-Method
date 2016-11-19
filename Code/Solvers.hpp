#ifndef SOLVERS_HPP
#define SOLVERS_HPP

#include<Eigen/Dense>
#include<tuple>


//Linear system for backward and forward substitution
Eigen::MatrixXd linearsystem_backward_subst(const Eigen::MatrixXd &U, const Eigen::MatrixXd &b);
Eigen::MatrixXd linearsystem_forward_subst(const Eigen::MatrixXd &U, const Eigen::MatrixXd &b);

//Linear Solver using LU decomposition with row pivoting
Eigen::VectorXd linear_solve_lu_row_pivoting(Eigen::MatrixXd &A, Eigen::VectorXd &b);

//Linear Solver using Cholesky decomposition
Eigen::VectorXd linear_solve_cholesky(Eigen::MatrixXd &A, Eigen::VectorXd &b);

//weekly percentage return(decreasing time order)
Eigen::MatrixXd weekly_percentage_return_d(Eigen::MatrixXd &A);

//weekly log return(decreasing time order)
Eigen::MatrixXd weekly_log_return_d(Eigen::MatrixXd &A);

//Covariance
Eigen::MatrixXd cov(Eigen::MatrixXd &A);

//Linear Squares Implementation
Eigen::MatrixXd linear_squares(Eigen::MatrixXd &A, Eigen::VectorXd &b);
#endif





