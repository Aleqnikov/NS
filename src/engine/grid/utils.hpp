#pragma once
#include "engine.hpp"


namespace grid_utils
{

    double get_theta(double max_density, std::size_t j);
    double get_radius(double q, double d_1, double R_s, std::size_t i);

    std::pair<double, double> get_progression_parameters(double v, double nu, double r, double r_border, std::size_t N);

}