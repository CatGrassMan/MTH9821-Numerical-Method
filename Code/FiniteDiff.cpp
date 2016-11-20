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
	double result = exp(tau - 2.0);
	return result;
}

double g_right(double tau)
{
	double result = exp(tau + 2.0);
	return result;
}

//Forward Eular(explicit method)
Eigen::MatrixXd forward_eular(int N, int M, double alpha,double x_left,double x_right,double t_final)
{
	Eigen::MatrixXd A = Eigen::MatrixXd::Zero(N - 1, N - 1);
	Eigen::VectorXd U = Eigen::VectorXd::Zero(N - 1);
	Eigen::VectorXd b = Eigen::VectorXd::Zero(N - 1);
	Eigen::MatrixXd valueAtNodes = Eigen::MatrixXd::Zero(M + 1, N + 1);
	double dx = (x_right - x_left) / N;
	double dt = t_final / M;
	//create matrix A
	for (int i = 0; i < N - 1; ++i)
	{
		A(i, i) = 1 - 2 * alpha;
	}
	for (int i = 0; i < N - 2; ++i)
	{
		A(i, i + 1) = alpha;
	}
	for (int i = 1; i < N - 1; ++i)
	{
		A(i, i - 1) = alpha;
	}
	//vector U
	for (int i = 0;i < N - 1; ++i)
	{
		U(i) = f(x_left + (i + 1)*dx);
	}
	//matrix valueAtNodes
	valueAtNodes(0, 0) = g_left(0);
	valueAtNodes(0, N) = g_right(0);
	for (int i = 1; i < N; ++i)
	{
		valueAtNodes(0, i) = U(i - 1);
	}
	for (int m = 0; m < M; ++m)
	{
		b(0) = alpha * g_left(m*dt);
		b(N - 2) = alpha * g_right(m*dt);
		U = A*U + b;
		valueAtNodes(m + 1, 0) = g_left((m + 1)*dt);
		valueAtNodes(m + 1, N) = g_right((m + 1)*dt);
		for (int i = 1; i < N; ++i)
		{
			valueAtNodes(m + 1, i) = U(i - 1);
		}
	}
	return valueAtNodes;
}

//Backward Eular(implicit method)
Eigen::MatrixXd backward_eular_lu(int N, int M, double alpha, double x_left, double x_right, double t_final)
{
	Eigen::MatrixXd A = Eigen::MatrixXd::Zero(N - 1, N - 1);
	Eigen::VectorXd U = Eigen::VectorXd::Zero(N - 1);
	Eigen::VectorXd b = Eigen::VectorXd::Zero(N - 1);
	Eigen::MatrixXd valueAtNodes = Eigen::MatrixXd::Zero(M + 1, N + 1);
	double dx = (x_right - x_left) / N;
	double dt = t_final / M;
	//create matrix A
	for (int i = 0; i < N - 1; ++i)
	{
		A(i, i) = 1 + 2 * alpha;
	}
	for (int i = 0; i < N - 2; ++i)
	{
		A(i, i + 1) = -alpha;
	}
	for (int i = 1; i < N - 1; ++i)
	{
		A(i, i - 1) = -alpha;
	}
	//vector U
	for (int i = 0; i < N - 1; ++i)
	{
		U(i) = f(x_left + (i + 1)*dx);
	}
	//matrix valueAtNodes
	valueAtNodes(0, 0) = g_left(0);
	valueAtNodes(0, N) = g_right(0);
	for (int i = 1; i < N; ++i)
	{
		valueAtNodes(0, i) = U(i - 1);
	}
	for (int m = 1; m <= M; ++m)
	{
		b(0) = alpha * g_left(m*dt);
		b(N - 2) = alpha * g_right(m*dt);
		Eigen::VectorXd U_sum = Eigen::VectorXd::Zero(N - 1);
		U_sum = U + b;
		U = linear_solve_lu_row_pivoting(A, U_sum);
		valueAtNodes(m, 0) = g_left(m*dt);
		valueAtNodes(m, N) = g_right(m*dt);
		for (int i = 1; i < N; ++i)
		{
			valueAtNodes(m, i) = U(i - 1);
		}
	}
	return valueAtNodes;

}
Eigen::MatrixXd backward_eular_sor(int N, int M, double alpha, double x_left, double x_right, double t_final, double w, double tol)
{
	Eigen::MatrixXd A = Eigen::MatrixXd::Zero(N - 1, N - 1);
	Eigen::VectorXd U = Eigen::VectorXd::Zero(N - 1);
	Eigen::VectorXd b = Eigen::VectorXd::Zero(N - 1);
	Eigen::MatrixXd valueAtNodes = Eigen::MatrixXd::Zero(M + 1, N + 1);
	double dx = (x_right - x_left) / N;
	double dt = t_final / M;
	//create matrix A
	for (int i = 0; i < N - 1; ++i)
	{
		A(i, i) = 1 + 2 * alpha;
	}
	for (int i = 0; i < N - 2; ++i)
	{
		A(i, i + 1) = -alpha;
	}
	for (int i = 1; i < N - 1; ++i)
	{
		A(i, i - 1) = -alpha;
	}
	//vector U
	for (int i = 0; i < N - 1; ++i)
	{
		U(i) = f(x_left + (i + 1)*dx);
	}
	//matrix valueAtNodes
	valueAtNodes(0, 0) = g_left(0);
	valueAtNodes(0, N) = g_right(0);
	for (int i = 1; i < N; ++i)
	{
		valueAtNodes(0, i) = U(i - 1);
	}
	for (int m = 1; m <= M; ++m)
	{
		b(0) = alpha * g_left(m*dt);
		b(N - 2) = alpha * g_right(m*dt);
		Eigen::VectorXd U_sum = Eigen::VectorXd::Zero(N - 1);
		U_sum = U + b;
		U = SOR_Iteration(A, U_sum, U, tol, w);
		valueAtNodes(m, 0) = g_left(m*dt);
		valueAtNodes(m, N) = g_right(m*dt);
		for (int i = 1; i < N; ++i)
		{
			valueAtNodes(m, i) = U(i - 1);
		}
	}
	return valueAtNodes;
}

//Crank Nicolson
Eigen::MatrixXd crank_nicolson_lu(int N, int M, double alpha, double x_left, double x_right, double t_final)
{
	Eigen::MatrixXd A = Eigen::MatrixXd::Zero(N - 1, N - 1);
	Eigen::MatrixXd B = Eigen::MatrixXd::Zero(N - 1, N - 1);
	Eigen::VectorXd U = Eigen::VectorXd::Zero(N - 1);
	Eigen::VectorXd b = Eigen::VectorXd::Zero(N - 1);
	Eigen::MatrixXd valueAtNodes = Eigen::MatrixXd::Zero(M + 1, N + 1);
	double dx = (x_right - x_left) / N;
	double dt = t_final / M;
	valueAtNodes(0, 0) = g_left(0);
	valueAtNodes(0, N) = g_right(0);
	//vector U
	for (int i = 0; i < N - 1; ++i)
	{
		U(i) = f(x_left + (i + 1)*dx);
	}
	for (int i = 1; i < N; ++i)
	{
		valueAtNodes(0, i) = U(i - 1);
	}
	//matrix A
	for (int i = 0; i < N - 1; ++i)
	{
		A(i, i) = 1.0 + alpha;
	}
	for (int i = 0; i < N - 2; ++i)
	{
		A(i, i + 1) = -alpha / 2.0;
	}
	for (int i = 1; i < N - 1; ++i)
	{
		A(i, i - 1) = -alpha / 2.0;
	}
	//matrix B
	for (int i = 0; i < N - 1; ++i)
	{
		B(i, i) = 1.0 - alpha;
	}
	for (int i = 0; i < N - 2; ++i)
	{
		B(i, i + 1) = alpha / 2.0;
	}
	for (int i = 1; i < N - 1; ++i)
	{
		B(i, i - 1) = alpha / 2.0;
	}
	//vector b
	Eigen::VectorXd U_sum1 = Eigen::VectorXd::Zero(N - 1);
	for (int m = 1; m <= M; ++m)
	{
		b(0) = alpha / 2.0 * (g_left(m*dt) + g_left((m -1)*dt));
		b(N - 2) = alpha / 2.0*(g_right(m*dt) + g_right((m - 1)*dt));
		U_sum1 = B*U + b;
		U = linear_solve_lu_row_pivoting(A, U_sum1);
		valueAtNodes(m, 0) = g_left(m*dt);
		valueAtNodes(m, N) = g_right(m*dt);
		for (int i = 1; i < N; ++i)
		{
			valueAtNodes(m, i) = U(i - 1);
		}
	}
	return valueAtNodes;
}
Eigen::MatrixXd crank_nicolson_sor(int N, int M, double alpha, double x_left, double x_right, double t_final, double w, double tol)
{
	Eigen::MatrixXd A = Eigen::MatrixXd::Zero(N - 1, N - 1);
	Eigen::MatrixXd B = Eigen::MatrixXd::Zero(N - 1, N - 1);
	Eigen::VectorXd U = Eigen::VectorXd::Zero(N - 1);
	Eigen::VectorXd b = Eigen::VectorXd::Zero(N - 1);
	Eigen::MatrixXd valueAtNodes = Eigen::MatrixXd::Zero(M + 1, N + 1);
	double dx = (x_right - x_left) / N;
	double dt = t_final / M;
	valueAtNodes(0, 0) = g_left(0);
	valueAtNodes(0, N) = g_right(0);
	//vector U
	for (int i = 0; i < N - 1; ++i)
	{
		U(i) = f(x_left + (i + 1)*dx);
	}
	for (int i = 1; i < N; ++i)
	{
		valueAtNodes(0, i) = U(i - 1);
	}
	//matrix A
	for (int i = 0; i < N - 1; ++i)
	{
		A(i, i) = 1.0 + alpha;
	}
	for (int i = 0; i < N - 2; ++i)
	{
		A(i, i + 1) = -alpha / 2.0;
	}
	for (int i = 1; i < N - 1; ++i)
	{
		A(i, i - 1) = -alpha / 2.0;
	}
	//matrix B
	for (int i = 0; i < N - 1; ++i)
	{
		B(i, i) = 1.0 - alpha;
	}
	for (int i = 0; i < N - 2; ++i)
	{
		B(i, i + 1) = alpha / 2.0;
	}
	for (int i = 1; i < N - 1; ++i)
	{
		B(i, i - 1) = alpha / 2.0;
	}
	//vector b
	Eigen::VectorXd U_sum1 = Eigen::VectorXd::Zero(N - 1);
	for (int m = 1; m <= M; ++m)
	{
		b(0) = alpha / 2.0 * (g_left(m*dt) + g_left((m - 1)*dt));
		b(N - 2) = alpha / 2.0*(g_right(m*dt) + g_right((m - 1)*dt));
		U_sum1 = B*U + b;
		U = SOR_Iteration(A, U_sum1, U, tol, w);
		valueAtNodes(m, 0) = g_left(m*dt);
		valueAtNodes(m, N) = g_right(m*dt);
		for (int i = 1; i < N; ++i)
		{
			valueAtNodes(m, i) = U(i - 1);
		}
	}
	return valueAtNodes;
}