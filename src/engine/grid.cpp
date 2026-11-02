#include "grid.hpp"

#include <iostream>

namespace
{

    double  get_theta(double max_density, std::size_t j)
    {
        return std::pi * j / max_density;
    }

    double get_radius(double q, double d_1, double R_s, std::size_t j)
    {
        return R_s + d_1 * (std::pow(q, i) - 1) / (q - 1);
    }

}


Grid::Grid(ModelConfig config)
{
    double height = config.radius * 20;

    int N_R_max = 1000, N_R_min = 100;
    int N_theta_max = 1500, N_theta_min = 150;

    int N_R_current  = N_R_min + config.density * (N_R_max - N_R_min);
    int N_theta_current = N_theta_min + config.density * (N_theta_max - N_theta_min);

    grid_.resize(N_R_current + 2, N_theta_current + 2);


    // Параметры самой сетки
    // коэфициент роста геометрической прогрессии.
    const double q = 1.0002;

    double d_1 = 19 * config.radius * (q  - 1) / (std::pow(q, N_R_current - 1) - 1);
    for (std::size_t i = 1; i < N_R_current + 1; ++i)
    {
        for (std::size_t j = 1; j < N_theta_current + 1; ++j)
        {
            double theta_last = get_theta(N_theta_current, j - 1);
            double theta_current = get_theta(N_theta_current, j);

            double R_last = get_radius(q, d_1, config.radius, i - 1);
            double R_current = get_radius(q, d_1, config.radius, i);

            grid_[i][j] = Cell(theta_current, theta_last, R_current, R_last, config.radius);

        }
    }

    this->init_cells();
    this->init_boundary_conditions(ModelConfig config);
}


std::pair<int, int> Grid::get_grid_dimensions()
{
    return {h_g, w_g};
}



CellGeometry::CellGeometry(double R_w, double R_e, double theta_s, double theta_n, double R_s)
{
	R_P = (R_w + R_e) / 2;
	theta_P = (theta_s  + theta_n) / 2;
	r_P  = R_P * std::sin(theta_P);
	d_P = R_P - R_s;
	delta_theta = theta_n - theta_s;
	delta_c = std::cos(theta_s) - std::cos(theta_n);
	V_P = (R_e * R_e * R_e - R_w * R_w * R_w) / 3  * delta_c;
	sigma_P  = (R_e * R_e - R_w * R_w ) / 2  * delta_theta;
}

EdgeGeometry::EdgeGeometry(CellGeometry gm)
{
	A_e = gm.R_e * gm.R_e * gm.delta_c;
	A_w = gm.R_w * gm.R_w *  gm.delta_c;
	A_n = std::sin(gm.theta_n) * (gm.R_e * gm.R_e - gm.R_w * gm.R_w )/ 2;
	A_s = std:sin(gm.theta_s) * (gm.R_e * gm.R_e - gm.R_w * gm.R_w )/ 2;

	S_z_e = gm.R_e * gm.R_e / 2 * (std::sin(gm.theta_n) * std::sin(gm.theta_n) - std::sin(gm.theta_s) * std::sin(gm.theta_s));
	S_r_e = gm.R_e * gm.R_e / 2  * (gm.delta_theta - (std::sin(2 * gm.theta_n)  - std::sin(2 * gm.theta_s)) / 2);

	S_z_w = gm.R_w * gm.R_w / 2 * (std::sin(gm.theta_n) * std::sin(gm.theta_n) - std::sin(gm.theta_s) * std::sin(gm.theta_s));
	S_r_w = gm.R_w * gm.R_w / 2  * (gm.delta_theta - (std::sin(2 * gm.theta_n)  - std::sin(2 * gm.theta_s)) / 2);

	S_z_n = - gm.A_n * std::sin(gm.theta_n);
	S_r_n = gm.A_n * std::cos(gm.theta_n);

	S_z_s = - gm.A_s * std::sin(gm.theta_n);
	S_r_s = gm.A_s *  std::cos(gm.theta_n);

}
