#include "Decomposition.hpp"
#include "Iteration.hpp"
#include "Substitution.hpp"
#include "Solvers.hpp"
#include "Interpolation.hpp"
#include "Optionpricing.hpp"
#include "FiniteDiff.hpp"
#include<iostream>
#include<Eigen/Dense>
#include<iomanip>
#include<tuple>
#include<math.h>
#include<vector>

//Funtion f and Funtion g_left and g_right
double f(double x)
{
	double result = exp(x);
	return result;
}

double g_left(double tau)
{
	double result = exp(tau - 2);
	return result;
}

double g_right(double tau)
{
	double result = exp(tau + 2);
	return result;
}

//Forward Eular(explicit method)
Eigen::VectorXd forward_eular(int N, int M, double alpha,double x_left,double x_right,double t_final)
{
	Eigen::MatrixXd A = Eigen::MatrixXd::Zero(N - 1, N - 1);
	Eigen::VectorXd U = Eigen::VectorXd::Zero(N - 1);
	Eigen::VectorXd b = Eigen::VectorXd::Zero(N - 1);
	double dx = (x_right - x_left) / N;
	double dt = t_final / M;
	for (int n = 0; n < N - 1; ++n)
	{
		U(n) = f(x_left + (n + 1)*dx);
		A(n, n) = 1 - 2 * alpha;
	}
	for (int n = 1; n < N - 1; ++n)
	{
		A(n, n - 1) = alpha;
		A(n - 1, n) = alpha;
	}
	for (int m = 0; m < M; ++m)
	{
		b(0) = g_left(m*dt);
		b(N - 2) = g_right(m*dt);
		U = A*U + b;
	}
	return U;
}