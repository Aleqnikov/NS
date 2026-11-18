#pragma once

#include "iter_ls_solver.hpp"
#include "util/solver_constants.hpp"

//https://algowiki-project.org/ru/%D0%A3%D1%87%D0%B0%D1%81%D1%82%D0%BD%D0%B8%D0%BA:%D0%A1%D0%BE%D1%80%D0%BE%D0%BA%D0%B8%D0%BD_%D0%90%D0%BB%D0%B5%D0%BA%D1%81%D0%B0%D0%BD%D0%B4%D1%80/%D0%9C%D0%B5%D1%82%D0%BE%D0%B4_%D1%81%D0%BE%D0%BF%D1%80%D1%8F%D0%B6%D0%B5%D0%BD%D0%BD%D1%8B%D1%85_%D0%B3%D1%80%D0%B0%D0%B4%D0%B8%D0%B5%D0%BD%D1%82%D0%BE%D0%B2_(%D0%A0%D0%B5%D1%88%D0%B5%D0%BD%D0%B8%D0%B5_%D0%A1%D0%9B%D0%90%D0%A3)
class CGSolver : public IterLsSolver {
public:
    CGSolver(int maxIterCount = constants::MAX_ITERATION_COUNT) : 
        LsSolverI(),  
        IterLsSolver(maxIterCount) {}
    
    std::vector<double> solve(const SparseMatrix& A, const std::vector<double>& b) override;
    std::vector<double> solve(const SparseMatrix& A, const std::vector<double>& b, std::vector<double>& x) override;
protected:
    bool solveLS(const SparseMatrix& A, const std::vector<double>& b, std::vector<double>& x) override;

    struct PriCondLC0 {
        //границы строк, если идти по матрице сначала по первой строке от начала до конца, потом по второй ...
        std::vector<int> row_bounds;
        //колонки cols[0] хранит в какой колонке находиться значение val[0]
        std::vector<int> cols;
        std::vector<double> vals;
        
        void build(const SparseMatrix& A);
        void precondition(const std::vector<double>& r, std::vector<double>& z) const;
    };
};
