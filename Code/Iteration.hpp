#ifndef ITERATION_HPP
#define ITERATION_HPP

#include<Eigen/Dense>
#include<tuple>

//Norm Function
double norm(Eigen::VectorXd &r);

//Jacobi iteration
Eigen::VectorXd jacobi_iteration(const Eigen::MatrixXd &A, const Eigen::VectorXd &b, Eigen::VectorXd &x0, double tol);

//Gauss Siedel iteration
Eigen::VectorXd GaussSiedel_Iteration(const Eigen::MatrixXd &A, const Eigen::VectorXd &b, Eigen::VectorXd &x0, double tol);

//SOR iteration
Eigen::VectorXd SOR_Iteration(const Eigen::MatrixXd &A, const Eigen::VectorXd &b, Eigen::VectorXd &x0, double tol, double w);

//SOR iteration return type is iteration
int SOR_Iteration_num(const Eigen::MatrixXd &A, const Eigen::VectorXd &b, Eigen::VectorXd &x0, double tol, double w);

#endif // !ITERATION_HPP



