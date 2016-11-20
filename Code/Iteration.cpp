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


//Norm function
double norm(Eigen::VectorXd &r)
{
	double sum = r.transpose()*r;
	double norm_number = sqrt(sum);
	return norm_number;
}

//Jacobi iteration
Eigen::VectorXd jacobi_iteration(const Eigen::MatrixXd &A, const Eigen::VectorXd &b, Eigen::VectorXd &x0, double tol)
{
	Eigen::MatrixXd L = Eigen::MatrixXd::Zero(A.rows(), A.cols());
	Eigen::MatrixXd U = Eigen::MatrixXd::Zero(A.rows(), A.cols());
	Eigen::MatrixXd D_Inverse = Eigen::MatrixXd::Zero(A.rows(), A.cols());
	for (int j = 0; j < A.rows(); ++j)
	{
		for (int k = 0; k < A.cols(); ++k)
		{
			if (j == k)
				D_Inverse(j, k) = 1 / A(j, k);
			else if (j < k)
				U(j, k) = A(j, k);
			else
				L(j, k) = A(j, k);

		}
	}
	Eigen::VectorXd x = x0;
	Eigen::VectorXd b_new = D_Inverse*b;
	int ic = 0;

	// Consecutive approximation criterion
	/*
	x0 = Eigen::VectorXd::Zero(A.rows());
	Eigen::VectorXd diff = x - x0;
	while (norm(diff)>tol)
	{
		x = -D_Inverse*(L+ U)*x + b_new;
		ic = ic + 1;
		diff = x - x0;
		x0 = x;
	}
	*/
	
	

	//Residue-Based criterion
	
	Eigen::VectorXd r0=b-A*x0;
	Eigen::VectorXd r=r0;
	double stop_iter_resid=tol*norm(r0);
	while(norm(r)>stop_iter_resid)
	{
	x=-D_Inverse*(L*x+U*x)+b_new;
	if (ic < 3)
		std::cout << "The " << ic << " times approximation is:" << std::endl<<x << std::endl << std::endl;
	r=b-A*x;
	ic=ic+1;
	}
	
	std::cout << "The iteration count of Jacobi is:"<< " " <<ic << std::endl << std::endl;
	return x;
}

//Gauss_Siedel iteration
Eigen::VectorXd GaussSiedel_Iteration(const Eigen::MatrixXd &A, const Eigen::VectorXd &b, Eigen::VectorXd &x0, double tol)
{
	Eigen::MatrixXd L = Eigen::MatrixXd::Zero(A.rows(), A.cols());
	Eigen::MatrixXd U = Eigen::MatrixXd::Zero(A.rows(), A.cols());
	Eigen::MatrixXd D = Eigen::MatrixXd::Zero(A.rows(), A.cols());
	for (int j = 0; j < A.rows(); ++j)
		for (int k = 0; k < A.cols(); ++k)
		{
			if (j == k)
				D(j, k) = A(j, k);
			else if (j < k)
				U(j, k) = A(j, k);
			else
				L(j, k) = A(j, k);
		}

	int ic = 0;
	Eigen::VectorXd x = x0;
	Eigen::VectorXd b_new = forward_subst(L + D, b);
	Eigen::MatrixXd R = linearsystem_forward_subst(L + D, U);

	// Consecutive approximation criterion
	/*
	x0=Eigen::VectorXd::Zero(A.rows());
	Eigen::VectorXd diff = x - x0;
	while(norm(diff)>tol)
	{
	x=-R*x+b_new;
	ic=ic+1;
	diff = x - x0;
	x0 = x;
	}
	*/
	

	//Residue-Based criterion
	
	Eigen::VectorXd r_0 = b - A*x0;
	Eigen::VectorXd r = r_0;
	double stop_iter_resid = tol*norm(r_0);
	while (norm(r) > stop_iter_resid)
	{
		x = -R*x + b_new;
		if (ic < 3)
			std::cout << "The " << ic << " times approximation is:" << std::endl << x << std::endl << std::endl;
		r = b - A*x;
		ic = ic + 1;
	}
	

	std::cout << "The iteration count of GS is:" << " " <<ic << std::endl << std::endl;
	return x;
}

//SOR iteration
Eigen::VectorXd SOR_Iteration(const Eigen::MatrixXd &A, const Eigen::VectorXd &b, Eigen::VectorXd &x0, double tol, double w)
{
	Eigen::MatrixXd L = Eigen::MatrixXd::Zero(A.rows(), A.cols());
	Eigen::MatrixXd U = Eigen::MatrixXd::Zero(A.rows(), A.cols());
	Eigen::MatrixXd D = Eigen::MatrixXd::Zero(A.rows(), A.cols());
	for (int j = 0; j < A.rows(); ++j)
		for (int k = 0; k < A.cols(); ++k)
		{
			if (j == k)
				D(j, k) = A(j, k);
			else if (j < k)
				U(j, k) = A(j, k);
			else
				L(j, k) = A(j, k);
		}
	int ic = 0;
	Eigen::VectorXd x = x0;
	Eigen::VectorXd r_0 = b - A*x0;
	Eigen::VectorXd r = r_0;
	double stop_iter_resid = tol*norm(r_0);
	Eigen::VectorXd b_new = w*forward_subst(w*L + D, b);
	Eigen::MatrixXd R = linearsystem_forward_subst(w*L + D, (1 - w)*D - w*U);

	//Consecutive criterion
	
	x0 = Eigen::VectorXd::Zero(A.rows());
	Eigen::VectorXd diff = x - x0;
	while (norm(diff)>tol)
	{
		x = R*x + b_new;
		ic = ic + 1;
		diff = x - x0;
		x0 = x;
	}
	
	

	//Residue-Based criterion
	/*
	while(norm(r) > stop_iter_resid)
	{
		x = R*x + b_new;
		if (ic < 3)
			std::cout << "The " << ic << " times approximation is:" << std::endl << x << std::endl << std::endl << std::endl;
		r = b - A*x;
		ic = ic + 1;
	}
	*/
	

	//std::cout << "The iteration count of SOR is:" << " " <<ic << std::endl << std::endl;
	return x;

}

//SOR iteration return type is iteration
int SOR_Iteration_num(const Eigen::MatrixXd &A, const Eigen::VectorXd &b, Eigen::VectorXd &x0, double tol, double w)
{
	Eigen::MatrixXd L = Eigen::MatrixXd::Zero(A.rows(), A.cols());
	Eigen::MatrixXd U = Eigen::MatrixXd::Zero(A.rows(), A.cols());
	Eigen::MatrixXd D = Eigen::MatrixXd::Zero(A.rows(), A.cols());
	for (int j = 0; j < A.rows(); ++j)
		for (int k = 0; k < A.cols(); ++k)
		{
			if (j == k)
				D(j, k) = A(j, k);
			else if (j < k)
				U(j, k) = A(j, k);
			else
				L(j, k) = A(j, k);
		}
	int ic = 0;
	Eigen::VectorXd x = x0;
	Eigen::VectorXd r_0 = b - A*x0;
	Eigen::VectorXd r = r_0;
	double stop_iter_resid = tol*norm(r_0);
	Eigen::VectorXd b_new = w*forward_subst(w*L + D, b);
	Eigen::MatrixXd R = linearsystem_forward_subst(w*L + D, (1 - w)*D - w*U);

	//Residue-Based criterion

	while (norm(r) > stop_iter_resid)
	{
		x = R*x + b_new;
		r = b - A*x;
		ic = ic + 1;
	}

	return ic;
}