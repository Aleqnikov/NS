#pragma once
#include "ls_solver_interface.hpp"

class IterLsSolver : public virtual LsSolverI {
public:
    IterLsSolver(int maxIterCount) : 
        maxIterCount(maxIterCount) {}
protected:
    virtual bool solveLS(const SparseMatrix& A, const std::vector<double>& b, std::vector<double>& x) = 0;
    virtual std::vector<double> solve(const SparseMatrix& A, const std::vector<double>& b, std::vector<double>& x) = 0;

    std::vector<double> calcResidual(const SparseMatrix& A, const std::vector<double>& b, 
        const std::vector<double>& x);

    //так называемый backward error - почитать как работает, я пока так и не понял
    double calcBackwardError(const SparseMatrix& A, const std::vector<double>& b, const std::vector<double>& x);

    //resTol - максимально терпимая невязка. Вычисляется после получения матриц.
    double resTol;
    const int maxIterCount;
};