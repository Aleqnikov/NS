#pragma once
#include "iter_ls_solver.hpp"
#include "util/solver_constants.h"

class GsSolver : public IterLsSolver {
public:
    GsSolver(
        double absDxTol = constants::ABSOLUTE_DX_TOLERANCE, 
        double relDxTol = constants::RELATIVE_DX_TOLERANCE,
        int maxIterCount = constants::MAX_ITERATION_COUNT,
        int errCheckInter = constants::ERROR_CHECK_ITER_INTERVAL) : 
        LsSolverI(), 
        IterLsSolver(maxIterCount), 
        absDxTol(absDxTol), 
        relDxTol(relDxTol),
        errCheckInter(errCheckInter) {}

    matrix solve(const std::vector<matrix>& A, const matrix& b) override;

    bool solveLS(const matrix& A, const std::vector<double>& b, std::vector<double>& x);

    double calcBackwardError(const matrix& A, const std::vector<double>& b, std::vector<double>& x);

protected:

    const double absDxTol;
    const double relDxTol;
    //вычисляется непосредствено перед началом решения системы
    double errTol = 0.0;
    const int errCheckInter;
};