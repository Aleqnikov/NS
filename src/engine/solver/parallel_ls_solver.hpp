#pragma once

#include "ls_solver_interface.hpp"

class ParallelLsSolver : public virtual LsSolverI {
public:
    ParallelLsSolver(int threadCount) : threadCount(threadCount) {};
protected:
    const int threadCount;
};