#pragma once
#include "ls_solver_interface.hpp"

class IterLsSolver : public virtual LsSolverI {
public:
    IterLsSolver(int maxIterCount) : 
        residualTolerance(residualTolerance), maxIterCount(maxIterCount) {}

    std::vector<double> calcResidual(const matrix& A, const std::vector<double>& b, 
        const std::vector<double>& x
    ) {
        std::vector<double> r;
        r.reserve(A.size());
        long double Ax_i;
        for(int i = 0; i < A.size(); i++) {
            Ax_i = 0;
            for(int j = 0; j < A.size(); j++) 
                Ax_i += (long double) A[i][j] * (long double) x[j];
            r.push_back((long double) b[i] - Ax_i);
        }
        return r;
    }

protected:
    //residualTolerance - максимально терпимая невязка. Вычисляется после получения матриц.
    const double residualTolerance;
    const int maxIterCount;
};