#pragma once
#ifndef OPTIONPRICING_HPP
#define OPTIONPRICING_HPP

//max function
double max(double a, double b);

//Binomial Tree Valuation of European Options
double BinomialTree_European(double S0, double sigma, double q, double r, double K, double N);

//Average Binomial Tree Valuation of Europen Options
double AverageBinomialTree_European(double S0, double sigma, double q, double r,double K, double N);

//Binomial Tree Valuation of American Put Options
double BinomialTree_American_Put(double S0, double sigma, double q, double r, double K, double N);

//Average Binomial Tree Valuation of American Options
double AverageBinomialTree_American(double S0, double sigma, double q, double r,  double K, double N);

#endif
