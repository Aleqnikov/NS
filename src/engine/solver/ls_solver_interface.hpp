#pragma once
#include <vector>
#include "../sparse_matrix.hpp"

class LsSolverI{
public:
    LsSolverI() = default;
    virtual ~LsSolverI() = default;

    /** Решает СЛАУ вида Ax=b
     * As - n матриц m*m (вроде как 5 на 5, но не принципиально) A[n][m][m]
     * bs - n матриц 1*m (тоже самое что одна матрица n*m) b[n][m]
     * return Xs - решение. n матриц 1*m (тоже самое что одна матрица n*m) X[n][m]
     */
    virtual std::vector<double> solve(const SparseMatrix& A, const std::vector<double>& b) = 0;
};