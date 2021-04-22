#include "newton.hpp"

void newton(dynamical_system* ds)
{
    Eigen::VectorXcd vp(ds->xdim+1);
    vp << ds->x0, ds->params(ds->var_param);
    // vp(Eigen::seqN(0, ds->xdim)) = ds->x0;
    // exit(0);
    std::cout << "******************************************" << std::endl;
    std::cout << "init = " << vp << std::endl;
    std::cout << "******************************************" << std::endl;
    Eigen::VectorXcd vn(ds->xdim+1);
    Eigen::VectorXcd F(ds->xdim+1);
    Eigen::MatrixXcd J(ds->xdim+1, ds->xdim+1);
    Eigen::MatrixXd jac(ds->xdim, ds->xdim);
    Eigen::VectorXcd eigvals(ds->xdim);
    Eigen::MatrixXcd I = Eigen::MatrixXd::Identity(ds->xdim, ds->xdim);
    // Eigen::VectorXcd temp(ds->xdim+1);
    double norm;

    std::ofstream file("out");
    
    for(int i = 0; i < ds->inc_iter; i++){
        for(int j = 0; j < ds->max_iter; j++){
            F = func_newton(vp, ds);
            jac = dTldx(vp(Eigen::seqN(0, ds->xdim)).real(), ds);
            J = jac_newton(vp, ds);
            
            vp = J.colPivHouseholderQr().solve(-F) + vp;
            norm = (vp-vn).norm();
            if(norm < ds->eps && j != 0){
                std::cout << "******************************************" << std::endl;
                std::cout << i << " : converged (iter = " << j << ")" << std::endl << vp << std::endl ;
                std::cout << "param = " << ds->params[ds->inc_param] << std::endl;
                Eigen::ComplexEigenSolver<Eigen::MatrixXd> eigensolver(jac);
                eigvals = eigensolver.eigenvalues();
                std::cout << eigvals << std::endl;
                std::cout << "abs = " << std::abs(eigvals(0)) << ",  arg = " << std::arg(eigvals(0)) << std::endl;
                std::cout << "******************************************" << std::endl;
                vn = vp;
                break;
            }

            if((vp-vn).norm() > ds->explode){
                std::cout << "explode" << std::endl;
                std::cout << j << std::endl;
                file.close();
                exit(1);
            }    

            if(j==ds->max_iter-1){
                std::cout << "iter over" << std::endl;
                std::cout << j << std::endl;
                file.close();
                exit(1);
            }
            
            file << ds->params[ds->inc_param] << " " << ds->params[ds->var_param] << std::endl;

            vn = vp;
        }
        ds->params[ds->inc_param] += ds->delta_inc;
    }
    file.close();
}