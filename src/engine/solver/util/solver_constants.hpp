#pragma once
#include <limits>

namespace constants
{
    static const double RELATIVE_DX_TOLERANCE = 1e-8;
    static const double ABSOLUTE_DX_TOLERANCE = RELATIVE_DX_TOLERANCE * 1e-2;
    static const int MAX_ITERATION_COUNT = 500;
    static const int ERROR_CHECK_ITER_INTERVAL = 10;
    //вычисление минимально допустимого значения ошибки, связаной с невязкой
    //где n - размерность системы. С - варируемый коэффициент [10, 100]
    static double calcErrorTolerance(int n) {
        int C = 100;
        return std::numeric_limits<double>::epsilon() * n * C;
    }
}