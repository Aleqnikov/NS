#include "utils.hpp"

#include <cmath>




namespace grid_utils


{
    const double PI = 3.14159265358979323846;



    double get_theta(double max_density, std::size_t j)
    {
        return PI * j / max_density;
    }

    double get_radius(double q, double d_1, double R_s, std::size_t i)
    {

        return R_s + d_1 * (std::pow(q, i) - 1) / (q - 1);
    }

    /** Сетка строится как геометрическая прогрессия. Первая ячейка расположена около шара, последняя
     * на расстоянии, заданном. Шаг сетки определяется размером самой первой ячейки, а ее размер, как раз таки
     * зависит от параметров самой системы.
     */
    std::pair<double, double> get_progression_parameters(double v, double nu, double r, double r_border, std::size_t N)
    {
        double Re = 2 * v * r / nu;
        double C_f = 0.026 / std::pow(Re, 1.0 / 7.0);
        double u = v * std::sqrt(C_f / 2);
        double d = 1 * nu / v;

        auto f = [&](double q)
        {
            return d * (std::pow(q, N) - 1) / (q - 1) - (r_border - r);
        };

        double lo = 1;
        double hi = 2;

        for (std::size_t i = 1; i <= 100; ++i)
        {
            double mid = (lo + hi) / 2;

            if (f(mid) > 0)
                hi = mid;
            else
                lo = mid;
        }

        double q = (lo + hi) / 2;

        return {d, q};
    }
}
