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



void main()
{
	std::cout.setf(std::ios::fixed);
	std::cout.precision(9);
	double x_left = -2.0;
	double x_right = 2.0;
	double t_final = 1.0;
	double alpha = 0.5;
	int M = 8;
	int N = 8;

	//forward eular
	Eigen::MatrixXd result1(M+1,N+1);
	result1 = forward_eular(N, M, alpha, x_left, x_right, t_final);
	std::cout << "The result of the forward eular is: " << std::endl;
	std::cout<<result1<<std::endl;

	//backward eular lu
	Eigen::MatrixXd result2(M + 1, N + 1);
	result2 = backward_eular_lu(N, M, alpha, x_left, x_right, t_final);
	std::cout << "The result of the backward eular lu is: " << std::endl;
	std::cout << result2 << std::endl;

	//backward eular sor
	double w = 1.2;
	double tol = 10E-6;
	Eigen::MatrixXd result3(M + 1, N + 1);
	result3 = backward_eular_sor(N, M, alpha, x_left, x_right, t_final, w, tol);
	std::cout << "The result of the backward eular sor is: " << std::endl;
	std::cout << result3 << std::endl;

	//crank nicolson lu
	Eigen::MatrixXd result4(M + 1, N + 1);
	result4 = crank_nicolson_lu(N, M, alpha, x_left, x_right, t_final);
	std::cout << "The result of the crank nicolson lu is: " << std::endl;
	std::cout << result4 << std::endl;

	//crank nicolson sor
	Eigen::MatrixXd result5(M + 1, N + 1);
	result5 = crank_nicolson_sor(N, M, alpha, x_left, x_right, t_final,w,tol);
	std::cout << "The result of the crank nicolson sor is: " << std::endl;
	std::cout << result5 << std::endl;

}