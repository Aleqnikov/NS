#include "CS_solver.hpp"
#include <iostream>
#include "util/util.hpp"
#include <cmath>
#include <limits>

std::vector<double> CGSolver::solve(const SparseMatrix& A, const std::vector<double>& b) {
    if (A.size() == 0)
        throw std::logic_error("A matrix can't empty");
    
    std::vector<double> x(A.size(), 0);
    resTol = constants::calcErrorTolerance(A.size());

    if(!solveLS(A, b, x)) 
        std::cout << "The matrix did not converge, the final error is: " << calcBackwardError(A, b, x) << std::endl;
    return x;
}


std::vector<double> CGSolver::solve(const SparseMatrix& A, const std::vector<double>& b, std::vector<double>& x) {
    if (A.size() == 0)
        throw std::logic_error("A matrix can't empty");
    
    resTol = constants::calcErrorTolerance(A.size());

    if(!solveLS(A, b, x)) 
        std::cout << "The matrix did not converge, the final error is: " << calcBackwardError(A, b, x) << std::endl;
    return x;
}


bool CGSolver::solveLS(const SparseMatrix& A, const std::vector<double>& b, std::vector<double>& x) {
    if(A.size() != x.size())
        x.resize(A.size(), 0);
    std::vector<double> r(b.size());
    {
        std::vector<double> Ax = util::multiply_sparse(A, x);
        for(int i = 0; i < b.size(); i ++) {
            r[i] = b[i] - Ax[i];
        }
    }
    
    std::vector<double> p = r;
    long double alph;
    double beta;


    double r_old_norm = util::multiply_vector(r, r);
    double r0_norm_sqrt = std::sqrt(r_old_norm);
    for(int k = 0; k < maxIterCount; k++) {
        if(std::abs(r_old_norm) < std::numeric_limits<double>::epsilon())
            return true;
        std::vector<double> Ap = util::multiply_sparse(A, p);
        double pAp =  util::multiply_vector(p, Ap);
        if(pAp <= 0) //может возникнуть только когда матрица не SPD
            return false;
        alph = r_old_norm / pAp;
        for(int i = 0; i < A.size(); i++) {
            x[i] = x[i] + alph * p[i];
            r[i] = r[i] - alph * Ap[i];
        }
        //проверка сходимости
        double r_new_norm = util::multiply_vector(r, r);
        if(std::sqrt(r_new_norm) < r0_norm_sqrt * resTol)
            return true;

        beta = r_new_norm / r_old_norm;
        r_old_norm = r_new_norm;
        for(int i = 0; i < A.size(); i++)
            p[i] = r[i] + beta * p[i];
    }
    return false;

}