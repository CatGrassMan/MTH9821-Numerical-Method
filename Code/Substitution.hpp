#ifndef SUBSTITUTION_HPP
#define SUBSTITUTION_HPP

#include<Eigen/Dense>
#include<tuple>


//Forward Substitution
Eigen::VectorXd forward_subst(const Eigen::MatrixXd &L, const Eigen::VectorXd &b);

//Backward Substitution
Eigen::VectorXd backward_subst(const Eigen::MatrixXd &U, const Eigen::VectorXd &b);


#endif // !SUBSTITUTION_HPP
