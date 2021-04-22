#include "newton.hpp"

void newton(dynamical_system* ds)
{
    Eigen::VectorXd vp = ds->x0;
    Eigen::VectorXd vn(ds->xdim);
    Eigen::VectorXd F(ds->xdim);
    Eigen::MatrixXd J(ds->xdim, ds->xdim);
    Eigen::MatrixXd jac(ds->xdim, ds->xdim);
    Eigen::VectorXcd eigvals;
    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(ds->xdim, ds->xdim);
    
    for(int i = 0; i < ds->inc_iter; i++){
        for(int j = 0; j < ds->max_iter; j++){
            F = Tl(vp, ds) - vp;
            jac = dTldx(vp, ds);
            J = jac - I;
            vp = J.colPivHouseholderQr().solve(-F) + vp;
            
            if((vp-vn).norm() < ds->eps && j != 0){
                std::cout << "******************************************" << std::endl;
                std::cout << i << " : converged (iter = " << j << ")" << std::endl;
                std::cout << std::fixed << std::setprecision(16) << vp << std::endl;
                std::cout << "param = ";
                std::cout << ds->params[ds->inc_param] << std::endl;
                std::cout << std::fixed << std::setprecision(8);
                Eigen::ComplexEigenSolver<Eigen::MatrixXd> eigensolver(jac);
                eigvals = eigensolver.eigenvalues();
                std::cout << eigvals << std::endl;
                std::cout << "norm = " << std::abs(eigvals[0]) << ",  arg = " << std::arg(eigvals[0]) << std::endl;
                std::cout << "******************************************" << std::endl;
                vn = vp;
                break;
            }

            if((vp-vn).norm() > ds->explode){
                std::cout << "explode" << std::endl;
                exit(1);
            }    

            if(j==ds->max_iter-1){
                std::cout << "iter over" << std::endl;
                exit(1);
            }

            vn = vp;
        }
        ds->params[ds->inc_param] += ds->delta_inc;
    }
}