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
	double x_left = -2;
	double x_right = 2;
	double t_final = 1;
	double alpha = 0.5;
	int M = 8;
	int N = 4;
	Eigen::VectorXd result(N - 1);
	result = forward_eular(N, M, alpha, x_left, x_right, t_final);
	std::cout << result;

	
}