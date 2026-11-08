#pragma once
#include "ls_solver_interface.hpp"

class IterLsSolver : public virtual LsSolverI {
public:
    IterLsSolver(double tolerance, int maxIterCount) : 
        tolerance(tolerance), maxIterCount(maxIterCount) {}

    /**
     * Можно передавать только матрицы формата
     * A - n*n
     * b - n*1
     * x - n*1
     */
    std::vector<double> calcResidual(const matrix& A, const std::vector<double>& b, 
        const std::vector<double>& x
    ) {
        std::vector<double> r;
        r.reserve(A.size());
        double Ax_i;
        for(int i = 0; i < A.size(); i++) {
            Ax_i = 0;
            for(int j = 0; j < A.size(); j++) 
                Ax_i += A[i][j] * x[j];
            r.push_back(b[i] - Ax_i);
        }
        return r;
    }

protected:
    //tolerance - максимально терпимая невязка
    const double tolerance;
    const int maxIterCount;
};