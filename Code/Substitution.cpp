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


//Forward Substitution
Eigen::VectorXd forward_subst(const Eigen::MatrixXd &L, const Eigen::VectorXd &b)
{
	Eigen::VectorXd x = Eigen::VectorXd::Zero(L.rows());
	x(0) = b(0) / L(0, 0);
	for (int j = 1; j < L.rows(); ++j)
	{
		double sum = 0;
		for (int k = 0; k < j; ++k)
			sum = sum + L(j, k)*x(k);
		x(j) = (b(j) - sum) / L(j, j);
	}
	return x;
}

//Backward Substitution
Eigen::VectorXd backward_subst(const Eigen::MatrixXd &U, const Eigen::VectorXd &b)
{
	Eigen::VectorXd x = Eigen::VectorXd::Zero(U.rows());
	int n = U.rows();
	x(n - 1) = b(n - 1) / U(n - 1, n - 1);
	for (int j = n - 2; j >= 0; --j)
	{
		double sum = 0;
		for (int k = j + 1; k < n; ++k)
			sum = sum + U(j, k)*x(k);
		x(j) = (b(j) - sum) / U(j, j);
	}
	return x;
}

