#include "dynamical_system.hpp"

dynamical_system::dynamical_system(nlohmann::json json){
    xdim = json["fixed"].size();
    this->period = json["period"];
    this->inc_param = json["inc_param"];
    this->var_param = json["var_param"];
    this->delta_inc = json["delta_inc"];
    this->inc_iter = json["inc_iter"];
    this->max_iter = json["max_iter"];
    this->eps = json["eps"];
    this->explode = json["explode"];

    /* These json array should be casted to the STL container type*/
    std::vector<double> fixed_arr = json["fixed"];
    Eigen::Map<Eigen::VectorXd> x0(fixed_arr.data(), fixed_arr.size());
    this->x0 = x0;

    std::vector<double> params_arr = json["params"];
    Eigen::Map<Eigen::VectorXd> params(params_arr.data(), params_arr.size());
    this->params = params;

    this->xk = std::vector<Eigen::VectorXd>(this->period);    
    std::vector<Eigen::MatrixXd> dTkdx(this->period);
    this->dTkdx = dTkdx;
}