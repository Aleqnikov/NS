#include "solver.hpp"

#include <iostream>
#include <ostream>
#include <cmath>


std::vector<double> GS (const std::vector<std::vector<double>>& A, std::vector<double> &b, double tol)
{
    std::size_t n = A.size();

    std::vector<double> x(n);

    for (std::size_t e = 0; e < 100; ++e)
    {
        std::vector<double> x_new = x;

        for (std::size_t i = 0; i < n; ++i)
        {
            double s_1 = 0;

            for (std::size_t j = 0; j < i; j++)
                s_1 += A[i][j] * x_new[j];

            double s_2 = 0;
            for (std::size_t j = i + 1; j < n; j++)
                s_2 += A[i][j] * x[j];

            x_new[i] = (b[i] - s_1 - s_2) / A[i][i];
        }
        x = x_new;


        double res = 0.0;
        for (std::size_t i = 0; i < n; ++i)
        {
            double r = b[i];
            for (std::size_t j = 0; j < n; ++j)
                r -= A[i][j] * x[j];

            res = std::max(res, std::fabs(r));
        }

        if (res < tol)
        {
            break;
        }

    }

    return x;
}

