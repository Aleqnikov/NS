#pragma once
#include <vector>

using matrix = std::vector<std::vector<double>>;
class LsSolverI{
public:
    LsSolverI() = default;
    ~LsSolverI() = default;

    /** Решает СЛАУ вида Ax=b
     * A - n матриц m*m (вроде как 5 на 5, но не принципиально) A[n][m][m]
     * b - n матриц 1*m (тоже самое что одна матрица n*m) b[n][m]
     * return X - решение. n матриц 1*m (тоже самое что одна матрица n*m) X[n][m]
     */
    virtual matrix solve(const std::vector<matrix>& A, const matrix& b) const = 0;
};