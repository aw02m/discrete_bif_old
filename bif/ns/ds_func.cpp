#include "ds_func.hpp"

// The operator[] is also overloaded for index-based access in vectors,
// but keep in mind that C++ doesn't allow operator[] to take more than one argument.
// We restrict operator[] to vectors, because an awkwardness in the C++ language
// would make matrix[i,j] compile to the same thing as matrix[j].

// Don't use std::complex<double> but Eigen::dcomplex.
// Eigen::dcomplex is just a wrapper of std::complex<doble>, however it is very flexible
// to substitute complex value to Eigen::Matrix classes.

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
    ret(1,0) = a * m * std::pow(x(0),m-1) * (x(0)*x(0)-b*b) +
               a * std::pow(x(0),m) * 2*x(0) + c;
    ret(1,1) = 0;
    return ret;
}

Eigen::VectorXd dTdlambda(Eigen::VectorXd x, dynamical_system* ds)
{
    Eigen::VectorXd ret(ds->xdim);
    double b = ds->params(1);
    unsigned int m = ds->params(4);

    switch(ds->var_param){
        case 0:
        ret(0) = 0;
        ret(1) = std::pow(x(0), m) * (x(0)*x(0) - b*b);
        break;
        case 3:
        ret(0) = x(1);
        ret(1) = 0;
        break;
    }

    return ret;
}

Eigen::MatrixXd dTdxdx(Eigen::VectorXd x, dynamical_system* ds, int k)
{
    Eigen::MatrixXd ret(ds->xdim, ds->xdim);
    double a, b, c, d;
    unsigned int m;
    a = ds->params(0);
    b = ds->params(1);
    c = ds->params(2);
    d = ds->params(3);
    m = ds->params(4);
    
    switch(k){
        case 0:
        ret(0,0) = 0;
        ret(0,1) = 0;
        ret(1,0) = a*m*(m-1)*std::pow(x(0), m-2)*(x(0)*x(0)-b*b) +
                   a*m*std::pow(x(0), m-1)*2*x(0) +
                   a*m*std::pow(x(0), m-1)*2*x(0) +
                   a*std::pow(x(0), m)*2;
        ret(1,1) = 0;
        break;
        
        case 1:
        ret(0,0) = 0;
        ret(0,1) = 0;
        ret(1,0) = 0;
        ret(1,1) = 0;
        break;
    }
    
    return ret;
}

Eigen::MatrixXd dTdxdlambda(Eigen::VectorXd x, dynamical_system* ds)
{
    Eigen::MatrixXd ret(ds->xdim, ds->xdim);
    double b = ds->params(1);
    unsigned int m = ds->params(4);

    switch(ds->var_param){
        case 0:
        ret(0,0) = 0;
        ret(0,1) = 0;
        ret(1,0) = m * std::pow(x(0), m-1) * (x(0)*x(0)-b*b) + 2 * std::pow(x(0), m+1);
        ret(1,1) = 0;
        break;
        case 3:
        ret(0,0) = 0;
        ret(0,1) = 1;
        ret(1,0) = 0;
        ret(1,1) = 0;
        break;
    }

    return ret;
}

//////////////////////
// DON'T EDIT BELOW //
//////////////////////

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

Eigen::VectorXd dTldlambda(Eigen::VectorXd x, dynamical_system* ds)
{
    Eigen::VectorXd ret(ds->xdim);

    for(int i = 0; i < ds->period; i++){
        ds->dTkdlambda[i] = dTdlambda(ds->xk[i], ds);
    }
    
    ret = ds->dTkdlambda[0];
    for(int i = 1; i < ds->period; i++){
        ret = ds->dTkdx[i] * ret + ds->dTkdlambda[i];
    }
    return ret;
}

Eigen::MatrixXd dTldxdx(Eigen::VectorXd x, dynamical_system* ds, int k)
{
    Eigen::MatrixXd ret(ds->xdim, ds->xdim);
    Eigen::MatrixXd A(ds->xdim, ds->xdim);
    Eigen::MatrixXd Ah(ds->xdim, ds->xdim);
    Eigen::VectorXd h_vec(ds->xdim);
    Eigen::VectorXd dummy(ds->xdim);
    double h = ds->dif_strip;

    for(int i = 0; i < ds->xdim; i++){
        if(i == k){
            h_vec(i) = h;
        }else{
            h_vec(i) = 0;
        }
    }
    
    dummy = Tl(x, ds);
    A = dTldx(x, ds);
    dummy = Tl(x+h_vec, ds);
    Ah = dTldx(x+h_vec, ds);

    return ret;
    // Eigen::MatrixXd ret(ds->xdim, ds->xdim);
    // Eigen::MatrixXd temp(ds->xdim, ds->xdim);
    // Eigen::MatrixXd O = Eigen::MatrixXd::Zero(ds->xdim, ds->xdim);
    // Eigen::MatrixXd I = Eigen::MatrixXd::Identity(ds->xdim, ds->xdim);
    // std::vector<Eigen::MatrixXd> prod_dTkdxdx(ds->period, I);
    // // std::cout << prod_dTkdx[1] << std::endl;

    // for(int i = 0; i < ds->period; i++){
    //     for(int j = i; j < ds->period; j++){
    //         if(i == j){
    //             // prod_dTkdx.push_back(prod_dTkdx[i] * dTdxdx(ds->xk[i], ds)[i]);
    //             prod_dTkdxdx[j] *= dTdxdx(ds->xk[i], ds, k);
    //         }else{
    //             prod_dTkdxdx[j] *= ds->dTkdx[j];
    //         }
    //     }
    // }
    
    // ret = O;
    // for(int i = 0; i < ds->period; i++){
    //     temp = I;
    //     for(int j = 0; j < ds->period; j++){
    //         if(i == j){
    //             // temp *= dTdxdx(ds->xk[ds->period-1-i], ds)[k];
    //             temp *= prod_dTkdxdx[i];
    //         }else{
    //             temp *= ds->dTkdx[ds->period-1-j];
    //         }
    //     }
    //     ret += temp;
    // }

    // return ret;
}

Eigen::MatrixXd dTldxdlambda(Eigen::VectorXd x, dynamical_system* ds)
{
    Eigen::MatrixXd ret(ds->xdim, ds->xdim);
    Eigen::MatrixXd A(ds->xdim, ds->xdim);
    Eigen::MatrixXd Ah(ds->xdim, ds->xdim);
    Eigen::VectorXd dummy(ds->xdim);
    double lambda = ds->params[ds->var_param];
    double h = ds->dif_strip;

    dummy = Tl(x, ds);
    A = dTldx(x, ds);
    ds->params[ds->var_param] += h;
    dummy = Tl(x, ds);
    Ah = dTldx(x, ds);
    ds->params[ds->var_param] = lambda;

    ret = (Ah - A)/h;

    return ret;
}


Eigen::dcomplex det_derivative(Eigen::MatrixXcd A, Eigen::MatrixXd dA, dynamical_system* ds)
{
    Eigen::MatrixXcd temp(ds->xdim, ds->xdim);
    Eigen::dcomplex ret(0,0);

    for(int i = 0; i < ds->xdim; i++){
        temp = A;
        temp.col(i) = dA.col(i).cast<Eigen::dcomplex>();
        ret += temp.determinant();
    }

    return ret;
}

Eigen::VectorXcd func_newton(Eigen::VectorXcd v, dynamical_system* ds)
{
    Eigen::VectorXd x = v(Eigen::seqN(0, ds->xdim)).real();
    ds->params[ds->var_param] = v(ds->xdim).real();
    Eigen::VectorXcd ret(ds->xdim+1);
    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(ds->xdim, ds->xdim);

    ret(Eigen::seqN(0, ds->xdim)) = Tl(x, ds) - x;

    Eigen::MatrixXd dtldx = dTldx(x, ds);
    Eigen::ComplexEigenSolver<Eigen::MatrixXd> eigensolver(dtldx);
    Eigen::VectorXcd eigvals = eigensolver.eigenvalues();
    ds->theta = std::arg(eigvals(0));
    Eigen::dcomplex mu(cos(ds->theta), sin(ds->theta));

    ret(ds->xdim) = (dtldx - mu*I).determinant();

    // std::cout << ret << std::endl;

    return ret;
}

Eigen::MatrixXcd jac_newton(Eigen::VectorXcd v, dynamical_system* ds)
{
    Eigen::VectorXd x = v(Eigen::seqN(0, ds->xdim)).real();
    ds->params[ds->var_param] = v(ds->xdim).real();
    Eigen::MatrixXcd ret(ds->xdim+1, ds->xdim+1);
    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(ds->xdim, ds->xdim);
    Eigen::MatrixXd dtldx = dTldx(x, ds);
    Eigen::dcomplex mu(cos(ds->theta), sin(ds->theta));
    Eigen::MatrixXcd A = dtldx - mu*I;

    ret(Eigen::seqN(0, ds->xdim), Eigen::seqN(0, ds->xdim)) = dtldx-I;
    ret(Eigen::seqN(0, ds->xdim), ds->xdim) = dTldlambda(x, ds);
    for(int i = 0; i < ds->xdim; i++){
        ret(ds->xdim, i) = det_derivative(A, dTdxdx(x, ds, i), ds);
    }
    ret(ds->xdim, ds->xdim) = det_derivative(A, dTldxdlambda(x, ds), ds);

    // std::cout << ret << std::endl;    

    return ret;
}