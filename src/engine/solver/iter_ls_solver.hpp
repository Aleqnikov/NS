#pragma once
#include "ls_solver_interface.hpp"

class IterLsSolver : public virtual LsSolverI {
public:
    IterLsSolver(double tolerance, int maxIterCount) : 
        tolerance(tolerance), maxIterCount(maxIterCount) {}

protected:
    const double tolerance;
    const int maxIterCount;
};