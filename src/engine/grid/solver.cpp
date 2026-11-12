#include "solver.hpp"

#include <iostream>
#include <ostream>
#include <cmath>




std::vector<double> GS(const SparseMatrix& A, const std::vector<double>& b, double tol)
{
    const std::size_t n = A.size();
    std::vector<double> x(n, 0.0);
    std::vector<double> x_new(n, 0.0);

    for (std::size_t e = 0; e < 7000; ++e)
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

        double diff_norm = 0.0;
        double x_norm    = 0.0;   //
        for (std::size_t i = 0; i < n; ++i)
        {
            diff_norm = std::max(diff_norm, std::fabs(x_new[i] - x[i]));
            x_norm    = std::max(x_norm,    std::fabs(x_new[i]));
        }
        const double rel_err = (x_norm > 0.0) ? (diff_norm / x_norm)
                                              : diff_norm;

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


        // остановка по относительной ошибке (и, для надёжности, по невязке)
        if (rel_err < tol ) break;
    }
    return x;
}