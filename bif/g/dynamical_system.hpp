#ifndef DYNAMICAL_SYSTEM_HPP_
#define DYNAMICAL_SYSTEM_HPP_

#include "sys_common.hpp"

class dynamical_system
{
    public:
    dynamical_system(nlohmann::json json);

    unsigned int xdim;

    Eigen::VectorXd x0;
    Eigen::VectorXd params;
    std::vector<Eigen::VectorXd> xk;

    std::vector<Eigen::MatrixXd> dTkdx;
    std::vector<Eigen::VectorXd> dTkdlambda;
    
    unsigned int period;
    unsigned int inc_param;
    unsigned int var_param;
    double delta_inc;
    unsigned int inc_iter;
    unsigned int max_iter;
    double eps;
    double explode;
    double dif_strip;
};

#endif