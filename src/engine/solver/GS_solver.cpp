#include "GS_solver.hpp"
#include <algorithm>
#include <stdexcept>
#include <iostream>

std::vector<double> GsSolver::solve(const SparseMatrix& A, const std::vector<double>& b) {
    if (A.size() == 0)
        throw std::logic_error("A matrix can't empty");
    
    std::vector<double> x(A.size(), 0);
    errTol = constants::calcErrorTolerance(A.size());
    
    if(!solveLS(A, b, x)) 
        std::cout << "The matrix did not converge, the final error is: " << calcBackwardError(A, b, x) << std::endl;
    return x;
}

std::vector<double> GsSolver::solve(const SparseMatrix& A, const std::vector<double>& b, std::vector<double>& x) {
    if (A.size() == 0)
        throw std::logic_error("A matrix can't empty");
    
    errTol = constants::calcErrorTolerance(A.size());
    
    if(!solveLS(A, b, x)) 
        std::cout << "The matrix did not converge, the final error is: " << calcBackwardError(A, b, x) << std::endl;
    return x;
}

bool GsSolver::solveLS(const SparseMatrix& A, const std::vector<double>& b, std::vector<double>& x) {
    x.resize(A.size());

    std::vector<double> A_diagonal(A.size(), 0);
    for(int i = 0; i < A.size(); i++) {
        for(int j = 0; j < A[i].indexes.size(); j++) {
            if(A[i].indexes[j] == i) {
                A_diagonal[i] = A[i].values[j];
                break;
            }
        }
    }

    //теперь проверяем на наличие значения на диагонали
    for(double val : A_diagonal) {
        if(val == 0)
            throw std::logic_error("zero or missing diagonal at row");
    }

    int k = 0;
    int dx_ind = -1;
    int checkResidualAcc = 0;
    double sum, dx, old_xi;
    while(k < maxIterCount) {
        dx = -1;
        for(int i = 0; i < A.size(); i++) {
            old_xi = x[i];
            sum = b[i];
            
            for(int j = 0; j < A[i].indexes.size(); j++) {
                if(A[i].indexes[j] == i) continue;
                sum -= A[i].values[j] * x[A[i].indexes[j]];
            }
            x[i] = sum / A_diagonal[i];

            if(dx < std::abs(old_xi - x[i])) {
                dx = std::abs(old_xi - x[i]);
                dx_ind = i;
            }
            
        }
        if (dx < absDxTol + std::abs(x[dx_ind]) * relDxTol)
            checkResidualAcc++;
        k++;

        if(checkResidualAcc >= errCheckInter) {
            if (calcBackwardError(A, b, x) < errTol) 
                return true;
            checkResidualAcc = 0;
        }
    }
    return false;
}


//так называемый backward error - почитать как работает, я пока так и не понял
double GsSolver::calcBackwardError(const SparseMatrix& A, const std::vector<double>& b, const std::vector<double>& x) {
    std::vector<double> residual = calcResidual(A, b, x);
    double A_norm = 0;
    for (const auto& row : A) {
        double s = 0.0;
        for (double v : row.values) s += std::abs(v);
        A_norm = std::max(A_norm, s);
    }
    double b_norm = 0, x_norm = 0, r_norm = 0;
    for (double v : b) b_norm = std::max(b_norm, std::abs(v));
    for (double v : x) x_norm = std::max(x_norm, std::abs(v));
    for (double v : residual) r_norm = std::max(r_norm, std::abs(v));
    double scale = A_norm * x_norm + b_norm;

    //только в случае если A = 0 и b = 0 (вроде как не возможно)
    if (scale == 0.0) return 0.0;
    return r_norm / scale;
}

