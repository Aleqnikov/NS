#include "grid.hpp"

#include <cmath>
#include <iostream>
#include <iomanip>

#include "utils.hpp"


namespace GU = grid_utils;


Grid::Grid(ModelConfig config)
{
    double height = config.radius * 20;

    int N_R_max = 500, N_R_min = 50;
    int N_theta_max = 2000, N_theta_min = 500;

    int N_R_current  = N_R_min + config.density * (N_R_max - N_R_min);
    int N_theta_current = N_theta_min + config.density * (N_theta_max - N_theta_min);

    grid_.resize(N_R_current + 2);
    for (auto& row : grid_)
        row.resize(N_theta_current + 2);

    auto res = GU::get_progression_parameters(config.speed, config.nu, config.radius, height, N_R_current);
    double d = res.first, q = res.second;

    for (std::size_t i = 1; i < N_R_current + 1; ++i)
    {
        for (std::size_t j = 1; j < N_theta_current + 1; ++j)
        {
            double theta_last = GU::get_theta(N_theta_current, j - 1);
            double theta_current = GU::get_theta(N_theta_current, j);

            double R_last = GU::get_radius(q, d, config.radius, i - 1);
            double R_current = GU::get_radius(q, d, config.radius, i);


            grid_[i][j] = Cell(theta_current, theta_last, R_current, R_last, config.radius);

        }
    }

}


std::pair<int, int> Grid::get_grid_dimensions()
{
    return {h_g, w_g};
}

Grid::Cell::Cell(double theta_0_, double theta_1_, double R_0_, double R_1_, double R_s)
    : theta_0(theta_0_), theta_1(theta_1_),
      R_0(R_0_), R_1(R_1_),
      geometry(R_1_, R_0_, theta_1_, theta_0_, R_s),
      edge_geometry(geometry)
{}

Grid::Cell::Cell(bool is_boundary_)
    : is_boundary(is_boundary_)
{}
