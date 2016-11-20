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


//Turn discount factor into zero rates
std::tuple<Eigen::VectorXd,Eigen::VectorXd> zero_rate(Eigen::VectorXd &t, Eigen::VectorXd &disc)
{
	Eigen::VectorXd time = Eigen::VectorXd::Zero(t.rows());
	for (int i = 0; i < t.rows(); ++i)
	{
		time(i) = t(i) / 12;
	}
	Eigen::VectorXd r = Eigen::VectorXd::Zero(disc.rows());
	for (int i = 0; i < disc.rows(); ++i)
	{
		r(i) = -log(disc(i)) / (t(i)/12);
	}
	std::tuple<Eigen::VectorXd, Eigen::VectorXd> result(time, r);
	return result;
}

//Efficient implementation of the natual cubic spline interpolation
//in the cubic spline, you have to add overnight rate into "x" and "v"
//x is the interpolation nodes and v is the interpolation value
std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd> efficient_cubic_spline(Eigen::VectorXd &x, Eigen::VectorXd &v)
{
	Eigen::VectorXd z = Eigen::VectorXd::Zero(x.rows() - 2);
	Eigen::MatrixXd M = Eigen::MatrixXd::Zero(x.rows() - 2, x.rows() - 2);
	for (int i = 1; i < x.rows() - 1; ++i)
	{
		z(i-1) = 6 * ((v(i + 1) - v(i)) / (x(i + 1) - x(i)) - (v(i) - v(i-1))/(x(i)-x(i-1)));
	}
	std::cout << "The vector z of the tridiagonal system is:" << std::endl << z << std::endl << std::endl;
	for (int i = 1; i < x.rows() - 1; ++i)
	{
		M(i-1, i-1) = 2 * (x(i + 1) - x(i-1));
	}
	for (int i = 1; i < x.rows() - 2; ++i)
	{
		M(i-1, i) = x(i + 1) - x(i);
	}
	for (int i = 2; i < x.rows() - 1; ++i)
	{
		M(i-1, i - 2) = x(i) - x(i-1);
	}
	std::cout << "The matrix M of the tridiagonal system is:" << std::endl << M << std::endl << std::endl;
	Eigen::VectorXd w = Eigen::VectorXd::Zero(x.rows() - 2);
	w = linear_solve_lu_row_pivoting(M, z);
	std::cout << "The linear solver result w is:" << std::endl << w << std::endl << std::endl;
	Eigen::VectorXd w_new = Eigen::VectorXd::Zero(x.rows());
	for (int i = 0; i < x.rows()-2; ++i)
	{
		w_new(i + 1) = w(i);
	}
	std::cout << "The n size linear solver result w_new is:" << std::endl << w_new << std::endl << std::endl;
	Eigen::VectorXd a = Eigen::VectorXd::Zero(x.rows() - 1);
	Eigen::VectorXd b = Eigen::VectorXd::Zero(x.rows() - 1);
	Eigen::VectorXd c = Eigen::VectorXd::Zero(x.rows() - 1);
	Eigen::VectorXd d = Eigen::VectorXd::Zero(x.rows() - 1);

	for (int i = 1; i < x.rows(); ++i)
	{
		c(i-1) = (w_new(i-1)*x(i) - w_new(i)*x(i-1)) / (2 * (x(i) - x(i-1)));
		d(i-1) = (w_new(i) - w_new(i-1)) / (6 * (x(i) - x(i-1)));
	}
	//std::cout << "The vector c is" << std::endl << c << std::endl << std::endl;
	//std::cout << "The vector d is" << std::endl << d << std::endl << std::endl;
	Eigen::VectorXd q = Eigen::VectorXd::Zero(x.rows() - 1);
	Eigen::VectorXd r= Eigen::VectorXd::Zero(x.rows() - 1);
	for (int i = 1; i < x.rows(); ++i)
	{
		q(i - 1) = v(i - 1) - c(i - 1)*x(i - 1)*x(i - 1) - d(i - 1)*x(i - 1)*x(i - 1)*x(i - 1);
		r(i-1) = v(i) - c(i - 1)*x(i)*x(i) - d(i - 1)*x(i)*x(i)*x(i);
	}
	//std::cout << "The vector q is" << std::endl << q << std::endl << std::endl;
	//std::cout << "The vector r is" << std::endl << r << std::endl << std::endl;
	for (int i = 1; i < x.rows(); ++i)
	{
		a(i - 1) = (q(i - 1)*x(i) - r(i-1)*x(i - 1)) / (x(i) - x(i - 1));
		b(i - 1) = (r(i-1) - q(i - 1)) / (x(i) - x(i - 1));
	}
	//std::cout << "The vector a is" << std::endl << a << std::endl << std::endl;
	//std::cout << "The vector b is" << std::endl << b << std::endl << std::endl;
	std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd> result(a,b,c,d);
	return result;
}


//Bond evaluation
double bond_value(Eigen::VectorXd &x, Eigen::VectorXd &v, int month,int frequency,double coupon_rate)
{
	std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd> result = efficient_cubic_spline(x, v);
	Eigen::VectorXd a = std::get<0>(result);
	Eigen::VectorXd b = std::get<1>(result);
	Eigen::VectorXd c = std::get<2>(result);
	Eigen::VectorXd d = std::get<3>(result);
	Eigen::VectorXd time = Eigen::VectorXd::Zero(month);
	Eigen::VectorXd z_rate = Eigen::VectorXd::Zero(month);
	for (int i = 0; i<month; ++i)
	{
		time(i) = (i+1)/12.0;
		for (int j = 0; j < x.rows() - 1; ++j)
		{
			if (time(i) <= x(j+1))
			{
				z_rate(i) = a(j) + b(j)*time(i) + c(j)*time(i)*time(i) + d(j)*time(i)*time(i)*time(i);
				break;
			}
		}
	}
	std::cout << "The time table for this bond is" << std::endl << time << std::endl << std::endl;
	std::cout << "The zero rate table for this bond is:" << std::endl << z_rate << std::endl << std::endl;
	
	//default principal = 100
	double principal = 100;
	double coupon = coupon_rate/(12/frequency)*principal;
	double sum = 0;
	int last = month - (int(month / frequency)) * frequency;
	for (int i=month-1; i>=(last-1); )
	{
		sum += coupon*exp(-time(i)*z_rate(i));
		i = i - frequency;
	}
	sum = sum + principal*exp(-time(month-1)*z_rate(month-1));
	return sum;
}
