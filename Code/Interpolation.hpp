#ifndef INTERPOLATION_HPP
#define INTERPOLATION_HPP

#include<Eigen/Eigen>

//Turn discount factor into zero rates
std::tuple<Eigen::VectorXd, Eigen::VectorXd> zero_rate(Eigen::VectorXd &t, Eigen::VectorXd &disc);

//Efficient implementation of the natual cubic spline interpolation
std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd> efficient_cubic_spline(Eigen::VectorXd &x, Eigen::VectorXd &v);

//Bond evaluation
double bond_value(Eigen::VectorXd &x, Eigen::VectorXd &v, int month, int frequency,double coupon_rate);

#endif
