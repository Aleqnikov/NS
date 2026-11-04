#include "grid.hpp"

#include <cmath>
#include <iostream>
#include <iomanip>

#include "utils.hpp"


namespace GU = grid_utils;


Grid::Grid(ModelConfig config) :
    config_(config)
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


            grid_[i][j] = Cell(theta_current, theta_last, R_current, R_last, config.radius, config);

        }
    }

    this->init_boundary_cond();

}


std::pair<int, int> Grid::get_grid_dimensions()
{
    return {h_g, w_g};
}

Grid::Cell::Cell(double theta_0_, double theta_1_, double R_0_, double R_1_, double R_s, ModelConfig config)
    : theta_0(theta_0_), theta_1(theta_1_),
      R_0(R_0_), R_1(R_1_),
      geometry(R_1_, R_0_, theta_1_, theta_0_, R_s),
      edge_geometry(geometry),
      fields(config.speed, config.nu)
{}

Grid::Cell::Cell(double theta_0_, double theta_1_, double R_0_, double R_1_, double R_s, ModelConfig config, bool is_boundary_)
    : Cell(theta_0_, theta_1_, R_0_, R_1_, R_s, config)
{
    is_boundary = is_boundary_;
}

Grid::Cell::Cell(bool is_boundary_): is_boundary(is_boundary_)
{
}

void Grid::init_boundary_cond()
{

    for (std::size_t i = 1; i < grid_.size() - 1; ++i)
    {
        std::size_t border = 0, upper = grid_[i].size() - 1;

        grid_[i][border] = Cell(true);
        grid_[i][border].fields.nut_P = grid_[i][border + 1].fields.nut_P;
        grid_[i][border].fields.p_P = grid_[i][border + 1].fields.p_P;
        grid_[i][border].fields.u_r_P = -grid_[i][border + 1].fields.u_r_P;
        grid_[i][border].fields.u_z_P = grid_[i][border + 1].fields.u_z_P;

        grid_[i][upper] = Cell(true);
        grid_[i][upper].fields.nut_P = grid_[i][upper - 1].fields.nut_P;
        grid_[i][upper].fields.p_P = grid_[i][upper - 1].fields.p_P;
        grid_[i][upper].fields.u_r_P = -grid_[i][upper - 1].fields.u_r_P;
        grid_[i][upper].fields.u_z_P = grid_[i][upper - 1].fields.u_z_P;
    }


    std::size_t left = 0, right = grid_.size() - 1;

    for (std::size_t j = 1; j < grid_[0].size() - 1; ++j)
    {
        grid_[left][j] = Cell(true);
        grid_[right][j] = Cell(true);

        // Стенка
        grid_[left][j].fields.nut_P = -grid_[left + 1][j].fields.nut_P;
        grid_[left][j].fields.p_P = grid_[left + 1][j].fields.p_P;
        grid_[left][j].fields.u_r_P = -grid_[left + 1][j].fields.u_r_P;
        grid_[left][j].fields.u_z_P = -grid_[left + 1][j].fields.u_z_P;


        if (grid_[right - 1][j].geometry.theta_P > M_PI / 2.0)
        {
            grid_[right][j].fields.nut_P = config_.nu * 6 - grid_[right - 1][j].fields.nut_P;
            grid_[right][j].fields.p_P = grid_[right - 1][j].fields.p_P;
            grid_[right][j].fields.u_r_P = -grid_[right - 1][j].fields.u_r_P;
            grid_[right][j].fields.u_z_P = config_.speed * 2 - grid_[right - 1][j].fields.u_z_P;
        }
        else
        {
            grid_[right][j].fields.nut_P = grid_[right - 1][j].fields.nut_P;
            grid_[right][j].fields.p_P = -grid_[right - 1][j].fields.p_P;
            grid_[right][j].fields.u_r_P = grid_[right - 1][j].fields.u_r_P;
            grid_[right][j].fields.u_z_P = grid_[right - 1][j].fields.u_z_P;
        }
    }
}

void Grid::update_boundary_cells()
{
    for (std::size_t i = 1; i < grid_.size() - 1; ++i)
    {
        std::size_t border = 0, upper = grid_[i].size() - 1;
        grid_[i][border].fields.nut_P = grid_[i][border + 1].fields.nut_P;
        grid_[i][border].fields.p_P = grid_[i][border + 1].fields.p_P;
        grid_[i][border].fields.u_r_P = -grid_[i][border + 1].fields.u_r_P;
        grid_[i][border].fields.u_z_P = grid_[i][border + 1].fields.u_z_P;

        grid_[i][upper].fields.nut_P = grid_[i][upper - 1].fields.nut_P;
        grid_[i][upper].fields.p_P = grid_[i][upper - 1].fields.p_P;
        grid_[i][upper].fields.u_r_P = -grid_[i][upper - 1].fields.u_r_P;
        grid_[i][upper].fields.u_z_P = grid_[i][upper - 1].fields.u_z_P;
    }


    std::size_t left = 0, right = grid_.size() - 1;

    for (std::size_t j = 1; j < grid_[0].size() - 1; ++j)
    {
        grid_[left][j].fields.nut_P = -grid_[left + 1][j].fields.nut_P;
        grid_[left][j].fields.p_P = grid_[left + 1][j].fields.p_P;
        grid_[left][j].fields.u_r_P = -grid_[left + 1][j].fields.u_r_P;
        grid_[left][j].fields.u_z_P = -grid_[left + 1][j].fields.u_z_P;


        if (grid_[right - 1][j].geometry.theta_P > M_PI / 2.0)
        {
            grid_[right][j].fields.nut_P = config_.nu * 6 - grid_[right - 1][j].fields.nut_P;
            grid_[right][j].fields.p_P = grid_[right - 1][j].fields.p_P;
            grid_[right][j].fields.u_r_P = -grid_[right - 1][j].fields.u_r_P;
            grid_[right][j].fields.u_z_P = config_.speed * 2 - grid_[right - 1][j].fields.u_z_P;
        }
        else
        {
            grid_[right][j].fields.nut_P = grid_[right - 1][j].fields.nut_P;
            grid_[right][j].fields.p_P = -grid_[right - 1][j].fields.p_P;
            grid_[right][j].fields.u_r_P = grid_[right - 1][j].fields.u_r_P;
            grid_[right][j].fields.u_z_P = grid_[right - 1][j].fields.u_z_P;
        }
    }
}
