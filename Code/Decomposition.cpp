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


//Matrix Multiplication
Eigen::MatrixXd matrix_multiply(const Eigen::MatrixXd &A, const Eigen::MatrixXd &B)
{
	Eigen::MatrixXd M = Eigen::ArrayXXd::Zero(A.rows(), B.cols());
	for (int j = 0; j<A.rows(); ++j)
		for (int k = 0; k<B.cols(); ++k)
		{
			// Inner product function: v1.dot(v2)
			Eigen::VectorXd v1 = A.block(j, 0, 1, A.cols()).transpose();
			Eigen::VectorXd v2 = B.block(0, k, B.rows(), 1);
			M(j, k) = v1.dot(v2);
		}
	return M;
}


//LU Decomposition without pivoting
std::tuple<Eigen::MatrixXd, Eigen::MatrixXd> lu_no_pivoting(Eigen::MatrixXd &A)
{
	// Initialize L with an identity matrix and U with a square matrixs with all the entries equal to 0
	Eigen::MatrixXd L = Eigen::MatrixXd::Identity(A.rows(), A.cols());
	Eigen::MatrixXd U = Eigen::ArrayXXd::Zero(A.rows(), A.cols());

	for (int j = 0; j<A.rows() - 1; ++j)
	{
		for (int k = j; k<A.rows(); ++k)
		{
			U(j, k) = A(j, k); // Calculate the jth row of U
			L(k, j) = A(k, j) / U(j, j); // Calculate the jth column of L
		}
		A.block(j + 1, j + 1, A.rows() - 1 - j, A.rows() - 1 - j) = A.block(j + 1, j + 1, A.rows() - 1 - j, A.rows() - 1 - j)
			- L.block(j + 1, j, A.rows() - 1 - j, 1)*U.block(j, j + 1, 1, A.rows() - 1 - j);
		// Update matrix A[j:n-1,j:n-1]
	}
	U(A.rows() - 1, A.cols() - 1) = A(A.rows() - 1, A.cols() - 1); // Calculate U[n-1,n-1]
	std::tuple<Eigen::MatrixXd, Eigen::MatrixXd> result(L, U);
	return result;
}

//LU Decomposition with row pivoting
std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd> lu_row_pivoting(Eigen::MatrixXd &A)
{
	// Initialize L and P with an identity matrix and U with a square matrixs with all the entries equal to 0
	Eigen::MatrixXd L = Eigen::MatrixXd::Identity(A.rows(), A.cols());
	Eigen::MatrixXd U = Eigen::MatrixXd::Zero(A.rows(), A.cols());
	Eigen::MatrixXd P = Eigen::MatrixXd::Identity(A.rows(), A.cols());

	for (int j = 0; j<A.rows() - 1; ++j)
	{
		// find the row index of the maximal absolute value within A[j:n-1,j] 
		int j_max = j;
		for (int k = j + 1; k<A.rows(); ++k)
		{
			if (abs(A(j_max, j))<abs(A(k, j)))
				j_max = k;
		}

		// Exchange j_max row and the jth row with column from j to n-1 of A
		Eigen::MatrixXd temp = A.block(j, j, 1, A.rows() - j);
		A.block(j, j, 1, A.rows() - j) = A.block(j_max, j, 1, A.rows() - j);
		A.block(j_max, j, 1, A.rows() - j) = temp;

		// Exchange j_max row and the jth row of P
		Eigen::MatrixXd temp1 = P.block(j, 0, 1, A.rows());
		P.block(j, 0, 1, A.rows()) = P.block(j_max, 0, 1, A.rows());
		P.block(j_max, 0, 1, A.rows()) = temp1;

		// Exchange the solved part of j_max row and j_th row of L
		if (j>0)
		{
			Eigen::MatrixXd temp3 = L.block(j, 0, 1, j);
			L.block(j, 0, 1, j) = L.block(j_max, 0, 1, j);
			L.block(j_max, 0, 1, j) = temp3;
		}

		for (int k = j; k<A.rows(); ++k)
		{
			U(j, k) = A(j, k); // Calculate the jth row of U
			L(k, j) = A(k, j) / U(j, j); // Calculate the jth column of L
		}
		A.block(j + 1, j + 1, A.rows() - 1 - j, A.rows() - 1 - j) = A.block(j + 1, j + 1, A.rows() - 1 - j, A.rows() - 1 - j)
			- L.block(j + 1, j, A.rows() - 1 - j, 1)*U.block(j, j + 1, 1, A.rows() - 1 - j);
		// Update matrix A[j:n-1,j:n-1]
	}
	U(A.rows() - 1, A.rows() - 1) = A(A.rows() - 1, A.rows() - 1); // Calculate U[n-1,n-1]
	std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd> result(P, L, U);
	return result;
}

//Cholesky
Eigen::MatrixXd cholesky(Eigen::MatrixXd &A)
{
	Eigen::MatrixXd U = Eigen::MatrixXd::Zero(A.rows(), A.cols());
	for(int i = 0; i < A.rows() - 1; ++i)
	{
		U(i, i) = sqrt(A(i, i));
		for (int k = i + 1; k < A.cols(); ++k)
		{
			U(i, k) = A(i, k) / U(i, i);
		}
		
		for (int j = i + 1; j < A.rows(); ++j)
		{
			for (int k = j; k < A.cols(); ++k)
			{
				A(j, k) = A(j, k) - U(i, j)*U(i, k);
			}
		}
	}
	U(A.rows() - 1, A.cols() - 1) = sqrt(A(A.rows() - 1, A.cols() - 1));
	return U;
}






