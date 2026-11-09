#pragma once

#include "GS_solver.hpp"
#include "parallel_ls_solver.hpp"
#include <atomic>

class ParallelGSSolver : public GsSolver, public ParallelLsSolver {
public:
    ParallelGSSolver(double tolerance = 0.01, double maxDx = 0.01, int maxIterCount = 100, int threadCount = 8):
        LsSolverI(), GsSolver(tolerance, maxDx, maxIterCount), ParallelLsSolver(threadCount) {}
    
    matrix solve(const std::vector<matrix>& As, const matrix& bs) override;

private:
    void worker(const std::vector<matrix>& As, const matrix& bs, matrix& Xs, std::atomic<int>& curInd, std::size_t threadId);
};