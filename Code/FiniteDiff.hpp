#pragma once
#ifndef FINITEDIFF_HPP
#define FINITEDIFF_HPP

#include<Eigen/Dense>

//Funtion g and Funtion f
double f(double x);
double g_left(double tau);
double g_right(double tau);

//Forward Eular(explicit method)
Eigen::MatrixXd forward_eular(int N, int M, double alpha,double x_left,double x_right,double t_final);

//Backward Eular(implicit method)
Eigen::MatrixXd backward_eular_lu(int N, int M, double alpha, double x_left, double x_right, double t_final);
Eigen::MatrixXd backward_eular_sor(int N, int M, double alpha, double x_left, double x_right, double t_final, double w, double tol);

//Crank Nicolson
Eigen::MatrixXd crank_nicolson_lu(int N, int M, double alpha, double x_left, double x_right, double t_final);
Eigen::MatrixXd crank_nicolson_sor(int N, int M, double alpha, double x_left, double x_right, double t_final, double w, double tol);



#endif

