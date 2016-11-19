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


//max function
double max(double a, double b)
{
	if (a > b)
	{
		return a;
	}
	else
	{
		return b;
	}
}

//Binomial Tree Valuation of European Put Options
double BinomialTree_European(double S0, double sigma, double q, double r, double K,double N)
{
	//We define MATURITY is 1
	double delta_t = 1.0 / N;
	double u = exp(sigma*sqrt(delta_t));
	double d = 1.0 / u;
	double Pu = (exp((r - q)*delta_t) - d) / (u - d);
	double Pd = 1.0 - Pu;
	std::vector<double> V(int(N) + 1);
	for (int i = 0; i < N + 1; ++i)
	{
		V[i] = max(K-S0*std::pow(u, N - i)*std::pow(d, i),0);
	}
	for (int j = int(N) - 1; j >= 0; --j)
	{
		for (int i = 0; i < j + 1; ++i)
		{
			V[i] = exp(-r*delta_t)*(Pu*V[i] + Pd*V[i + 1]);
		}
	}
	double V_final;
	V_final = V[0];
	return V_final;
}

//Average Binomial Tree Valuation of Europen Options
double AverageBinomialTree_European(double S0, double sigma, double q, double r, double K, double N)
{
	double V_bino_N = BinomialTree_European(S0, sigma, q, r, K, N);
	double V_bino_N_1 = BinomialTree_European(S0, sigma, q, r, K, N + 1);
	double V_average = 0.5*(V_bino_N + V_bino_N_1);
	return V_average;
}

//Binomial Tree Valuation of American Put Options
double BinomialTree_American_Put(double S0, double sigma, double q, double r, double K, double N)
{
	//We define MATURITY is 1
	double delta_t = 1.0 / N;
	double u = exp(sigma*sqrt(delta_t));
	double d = 1.0 / u;
	double Pu = (exp((r - q)*delta_t) - d) / (u - d);
	double Pd = 1.0 - Pu;
	std::vector<double> V(int(N) + 1);
	for (int i = 0; i < N + 1; ++i)
	{
		V[i] = max(K - S0*pow(u, N - i)*pow(d, i), 0);
	}
	for (int j = int(N) - 1; j >= 0; --j)
	{
		for (int i = 0; i < j + 1; ++i)
		{
			V[i] = max(exp(-r*delta_t)*(Pu*V[i] + Pd*V[i + 1]), K - S0*pow(u, j - i)*pow(d, i));
		}
	}
	double V_final;
	V_final = V[0];
	return V_final;
}

//Average Binomial Tree Valuation of American Options
double AverageBinomialTree_American(double S0, double sigma, double q, double r, double K, double N)
{
	double V_bino_N = BinomialTree_American_Put(S0, sigma, q, r, K, N);
	double V_bino_N_1 = BinomialTree_American_Put(S0, sigma, q, r, K, N + 1);
	double V_average = 0.5*(V_bino_N + V_bino_N_1);
	return V_average;
}