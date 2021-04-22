#ifndef DS_FUNC_HPP_
#define DS_FUNC_HPP_

#include "sys_common.hpp"

class dynamical_system;

Eigen::VectorXd T(Eigen::VectorXd x, dynamical_system* ds);
Eigen::VectorXd Tl(Eigen::VectorXd x, dynamical_system* ds);
Eigen::MatrixXd dTdx(Eigen::VectorXd x, dynamical_system* ds);
Eigen::MatrixXd dTldx(Eigen::VectorXd x, dynamical_system* ds);

#endif