#include "ds_func.hpp"

// The operator[] is also overloaded for index-based access in vectors,
// but keep in mind that C++ doesn't allow operator[] to take more than one argument.
// We restrict operator[] to vectors, because an awkwardness in the C++ language
// would make matrix[i,j] compile to the same thing as matrix[j].

Eigen::VectorXd T(Eigen::VectorXd x, dynamical_system* ds)
{
    Eigen::VectorXd ret(ds->xdim);
    double a, b, c, d;
    unsigned int m;
    a = ds->params(0);
    b = ds->params(1);
    c = ds->params(2);
    d = ds->params(3);
    m = ds->params(4);

    ret(0) = d * x(1);
    ret(1) = a * std::pow(x(0), m) * (x(0)*x(0) - b*b) + c*x(0);
    
    return ret;
}

Eigen::MatrixXd dTdx(Eigen::VectorXd x, dynamical_system* ds)
{
    Eigen::MatrixXd ret(ds->xdim, ds->xdim);
    double a, b, c, d;
    unsigned int m;
    a = ds->params(0);
    b = ds->params(1);
    c = ds->params(2);
    d = ds->params(3);
    m = ds->params(4);

    ret(0,0) = 0;
    ret(0,1) = d;
    ret(1,0) = a * m * std::pow(x(0),m-1) * (x(0)*x(0)-b*b) + a * std::pow(x(0),m) * 2*x(0) + c;
    ret(1,1) = 0;

    return ret;
}

Eigen::VectorXd Tl(Eigen::VectorXd x, dynamical_system* ds)
{
    ds->x0 = x;
    ds->xk[0] = ds->x0;
    
    for(int i = 1; i < ds->period; i++){
        ds->xk[i] = T(ds->xk[i-1], ds);
    }
    
    return T(ds->xk[ds->period-1], ds);
}


Eigen::MatrixXd dTldx(Eigen::VectorXd x, dynamical_system* ds)
{
    Eigen::MatrixXd ret = Eigen::MatrixXd::Identity(ds->xdim, ds->xdim);

    for(int i = 0; i < ds->period; i++){
        ds->dTkdx[i] = dTdx(ds->xk[i], ds);
    }

    for(int i = ds->period-1; i >= 0; i--){
        ret *= ds->dTkdx[i];
    }

    return ret;
}