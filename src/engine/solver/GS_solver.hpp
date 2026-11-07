#pragma once
#include "iter_ls_solver.hpp"

class GsSolver : public IterLsSolver {
public:
    GsSolver(double tolerance = 0.01, int maxIterCount = 100) : 
        LsSolverI(), IterLsSolver(tolerance, maxIterCount) {}

    matrix solve(const std::vector<matrix>& A, const matrix& b) const override; 
};