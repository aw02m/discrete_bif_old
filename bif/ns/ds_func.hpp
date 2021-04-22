#ifndef DS_FUNC_HPP_
#define DS_FUNC_HPP_

#include "sys_common.hpp"

class dynamical_system;

Eigen::VectorXd T(Eigen::VectorXd x, dynamical_system* ds);
Eigen::VectorXd Tl(Eigen::VectorXd x, dynamical_system* ds);
Eigen::MatrixXd dTdx(Eigen::VectorXd x, dynamical_system* ds);
Eigen::MatrixXd dTldx(Eigen::VectorXd x, dynamical_system* ds);
Eigen::VectorXd dTdlambda(Eigen::VectorXd x, dynamical_system* ds);
Eigen::VectorXd dTldlambda(Eigen::VectorXd x, dynamical_system* ds);
Eigen::MatrixXd dTdxdx(Eigen::VectorXd x, dynamical_system* ds, int k);
Eigen::MatrixXd dTldxdx(Eigen::VectorXd x, dynamical_system* ds, int k);
Eigen::MatrixXd dTdxdlambda(Eigen::VectorXd x, dynamical_system* ds);
Eigen::MatrixXd dTldxdlambda(Eigen::VectorXd x, dynamical_system* ds);
double det_derivative(Eigen::MatrixXd A, Eigen::MatrixXd dA, dynamical_system* ds);
Eigen::VectorXcd func_newton(Eigen::VectorXcd v, dynamical_system* ds);
Eigen::MatrixXcd jac_newton(Eigen::VectorXcd v, dynamical_system* ds);

#endif