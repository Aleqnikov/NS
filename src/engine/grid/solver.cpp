#include "solver.hpp"

#include <iostream>
#include <ostream>
#include <cmath>


std::vector<double> GS(const SparseMatrix& A, const std::vector<double>& b, double tol)
{
    const std::size_t n = A.size();
    std::vector<double> x(n, 0.0);
    std::vector<double> x_new(n, 0.0);

    for (std::size_t e = 0; e < 10000; ++e)
    {
        // прямой ход Гаусса-Зейделя
        for (std::size_t i = 0; i < n; ++i)
        {
            double s_1 = 0.0, s_2 = 0.0, diag = 0.0;
            const auto& row = A[i];
            for (std::size_t k = 0; k < row.nnz; ++k)
            {
                const std::size_t j   = row.entries[k].first;
                const double      val = row.entries[k].second;
                if (j < i)       s_1  += val * x_new[j];
                else if (j > i)  s_2  += val * x[j];
                else             diag  = val;
            }
            x_new[i] = (b[i] - s_1 - s_2) / diag;
        }
        x.swap(x_new);

        // невязка — тоже по nnz, а не по n
        double res = 0.0;
        for (std::size_t i = 0; i < n; ++i)
        {
            double r = b[i];
            const auto& row = A[i];
            for (std::size_t k = 0; k < row.nnz; ++k)
                r -= row.entries[k].second * x[row.entries[k].first];
            res = std::max(res, std::fabs(r));
        }
        if (res < 1e-8) break;
    }
    return x;
}
