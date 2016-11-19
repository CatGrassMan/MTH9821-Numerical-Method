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


//Linear system for backward substitution
Eigen::MatrixXd linearsystem_backward_subst(const Eigen::MatrixXd &U, const Eigen::MatrixXd &b)
{
	Eigen::MatrixXd X = Eigen::ArrayXXd::Zero(b.rows(), b.cols());
	for (int k = 0; k<X.cols(); ++k)
	{
		X.block(0, k, X.rows(), 1) = backward_subst(U, b.block(0, k, b.rows(), 1));
	}
	return X;

}

//Linear system for forward substitution
Eigen::MatrixXd linearsystem_forward_subst(const Eigen::MatrixXd &L, const Eigen::MatrixXd &b)
{
	Eigen::MatrixXd X = Eigen::ArrayXXd::Zero(b.rows(), b.cols());
	for (int k = 0; k<X.cols(); ++k)
	{
		X.block(0, k, X.rows(), 1) = forward_subst(L, b.block(0, k, b.rows(), 1));
	}
	return X;

}

//Linear Solver using LU decomposition with row pivoting
Eigen::VectorXd linear_solve_lu_row_pivoting(Eigen::MatrixXd &A, Eigen::VectorXd &b)
{
	std::tuple<Eigen::MatrixXd, Eigen::MatrixXd,Eigen::MatrixXd> result = lu_row_pivoting(A);
	Eigen::VectorXd b_new = std::get<0>(result)*b;
	Eigen::VectorXd y = forward_subst(std::get<1>(result), b_new);
	Eigen::VectorXd x = backward_subst(std::get<2>(result), y);
	return x;
}

//Linear Solver using Cholesky decomposition
Eigen::VectorXd linear_solve_cholesky(Eigen::MatrixXd &A, Eigen::VectorXd &b)
{
	Eigen::MatrixXd U = cholesky(A);
	Eigen::MatrixXd U_transpose = U.transpose();
	Eigen::VectorXd y = forward_subst(U_transpose, b);
	Eigen::VectorXd x = backward_subst(U, y);
	return x;
}

//weekly percentage return(decreasing time order)
Eigen::MatrixXd weekly_percentage_return_d(Eigen::MatrixXd &A)
{
	Eigen::MatrixXd weekly_return = Eigen::MatrixXd::Zero(A.rows() - 1, A.cols());
	for (int i = 0; i < A.rows()-1; ++i)
	{
		for (int j = 0; j < A.cols(); ++j)
		{
			weekly_return(i, j) = (A(i, j) - A(i+1, j)) / A(i+1, j);
		}
	}
	return weekly_return;
}

//weekly log return(decreasing time order)
Eigen::MatrixXd weekly_log_return_d(Eigen::MatrixXd &A)
{
	Eigen::MatrixXd weekly_return = Eigen::MatrixXd::Zero(A.rows() - 1, A.cols());
	for (int i = 0; i < A.rows() - 1; ++i)
	{
		for (int j = 0; j < A.cols(); ++j)
		{
			weekly_return(i, j) = log(A(i, j) / A(i + 1, j));
		}
	}
	return weekly_return;
}

//Covariance
Eigen::MatrixXd cov(Eigen::MatrixXd &A)
{
	Eigen::VectorXd mean = Eigen::VectorXd::Zero(A.cols());
	for (int i = 0; i < A.cols(); i++)
		mean(i) = A.col(i).mean();
	Eigen::VectorXd one = Eigen::VectorXd::Ones(A.rows());
	Eigen::MatrixXd Covariance = ((A - one*(mean.transpose())).transpose())*(A - one*(mean.transpose())) / (A.rows() - 1);
	return Covariance;
}

//Linear Squares Implementation
Eigen::MatrixXd linear_squares(Eigen::MatrixXd &A, Eigen::VectorXd &y)
{
	Eigen::MatrixXd A_tA = A.transpose()*A;
	Eigen::VectorXd A_ty = A.transpose()*y;
	Eigen::VectorXd x = linear_solve_cholesky(A_tA, A_ty);
	return x;
}