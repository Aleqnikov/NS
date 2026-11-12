#include "GS_solver.hpp"
#include <algorithm>
#include <stdexcept>

matrix GsSolver::solve(const std::vector<matrix>& A, const matrix& b) {
    matrix Xs(A.size());
    errTol = constants::calcErrorTolerance(A.size());
    
}

bool GsSolver::solveLS(const matrix& A, const std::vector<double>& b, std::vector<double>& x) {
    x.resize(A.size());
    for(int i = 0; i < A.size(); i++) {
        x[i] = b[i]/A[i][i];
    }

    int k = 0;
    int dx_ind = -1;
    int needToCheckResidualAcc = 0;
    double sum, dx, old_xi;
    while(k < maxIterCount) {
        dx = -1;
        for(int i = 0; i < A.size(); i++) {
            old_xi = x[i];
            sum = b[i];
            for(int j = 0; j < i; j++) {
                sum -= A[i][j] * x[j];
            }
            for(int j = i+1; j < A.size(); j++) {
                sum -= A[i][j] * x[j];
            }
            x[i] = sum / A[i][i];
            if(dx < std::abs(old_xi - x[i])) {
                dx = std::abs(old_xi - x[i]);
                dx_ind = i;
            }
            
        }
        if (dx < absDxTol + x[dx_ind] * relDxTol)
            needToCheckResidualAcc++;
        k++;

        if(needToCheckResidualAcc >= errCheckInter) {
            if (calcBackwardError(A, b, x) < errTol) 
                return true;
            needToCheckResidualAcc = 0;
        }
    }
    return false;
}


//так называемый backward error - почитать как работает, я пока так и не понял
double GsSolver::calcBackwardError(const matrix& A, const std::vector<double>& b, std::vector<double>& x) {
    std::vector<double> residual = calcResidual(A, b, x);
    double A_norm = 0;
    for (const auto& row : A) {
        double s = 0.0;
        for (double v : row) s += std::abs(v);
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

