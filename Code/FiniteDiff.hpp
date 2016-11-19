#pragma once
#ifndef FINITEDIFF_HPP
#define FINITEDIFF_HPP

#include<Eigen/Dense>

//Funtion g and Funtion f
double f(double x);
double g_left(double tau);
double g_right(double tau);

//Forward Eular(explicit method)
Eigen::VectorXd forward_eular(int N, int M, double alpha,double x_left,double x_right,double t_final);


#endif

