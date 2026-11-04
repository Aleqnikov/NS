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


void Grid::Cell::update_SA_parameters(Cell E, Cell W, Cell N, Cell S, ModelConfig config)
{


    sa_params.xi_P = fields.nut_P / config.nu;
    sa_params.f_v1_P = std::pow(sa_params.xi_P, 3) / (std::pow(sa_params.xi_P, 3) + std::pow(config.c_v_1, 3));
    sa_params.f_v2_P = 1 - sa_params.xi_P / (1 + sa_params.xi_P * sa_params.f_v1_P);
    sa_params.nu_t_P = fields.nut_P * sa_params.f_v1_P;
    sa_params.nu_eff_P = config.nu + sa_params.nu_t_P;
    sa_params.gamma_P = config.nu + fields.nut_P;

    double du_zd_R =  (E.fields.u_z_P - W.fields.u_z_P)/(E.geometry.R_P - W.geometry.R_P);
    double du_zd_theta = (N.fields.u_r_P - W.fields.u_r_P) / (N.geometry.theta_P - W.geometry.theta_P);


    double du_zdr = std::sin(geometry.theta_P) * du_zd_R +
        std::cos(geometry.theta_P)/geometry.R_P * du_zd_theta;


    double du_rd_R =  (E.fields.u_r_P - W.fields.u_r_P)/(E.geometry.R_P - W.geometry.R_P);
    double du_rd_theta = (N.fields.u_r_P - W.fields.u_r_P) / (N.geometry.theta_P - W.geometry.theta_P);

    double du_r_d_theta = std::cos(geometry.theta_P) * du_rd_R  - std::sin(geometry.theta_P) / geometry.R_P * du_rd_theta;

    sa_params.Omega_P = std::fabs(du_zdr - du_r_d_theta);

    sa_params.S_tilda_P = sa_params.Omega_P + fields.nut_P/(config.k * config.k * (geometry.d_P * geometry.d_P)) * sa_params.f_v2_P;

    sa_params.r_sa_P = std::min(10., fields.nut_P/(config.k * config.k * (geometry.d_P * geometry.d_P) * sa_params.S_tilda_P));
    sa_params.g_P = sa_params.r_sa_P + config.c_w_2 * (std::pow(sa_params.r_sa_P, 6) - sa_params.r_sa_P);
    sa_params.f_w_P = sa_params.g_P * std::pow(((1 + std::pow(config.c_w_3, 6))/(std::pow(sa_params.g_P, 6) + std::pow(config.c_w_3, 6))), 1.0/6.0);

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

        grid_[i][border].edge_geometry.F_e = 0;
        grid_[i][border].edge_geometry.F_s = 0;
        grid_[i][border].edge_geometry.F_n = 0;
        grid_[i][border].edge_geometry.F_w = 0;

        grid_[i][upper].edge_geometry.F_e = 0;
        grid_[i][upper].edge_geometry.F_s = 0;
        grid_[i][upper].edge_geometry.F_n = 0;
        grid_[i][upper].edge_geometry.F_w = 0;
    }


    std::size_t left = 0, right = grid_.size() - 1;

    for (std::size_t j = 1; j < grid_[0].size() - 1; ++j)
    {
        grid_[left][j] = Cell(true);
        grid_[right][j] = Cell(true);

        grid_[left][j].fields.nut_P = -grid_[left + 1][j].fields.nut_P;
        grid_[left][j].fields.p_P = grid_[left + 1][j].fields.p_P;
        grid_[left][j].fields.u_r_P = -grid_[left + 1][j].fields.u_r_P;
        grid_[left][j].fields.u_z_P = -grid_[left + 1][j].fields.u_z_P;


        grid_[left][j].edge_geometry.F_e = 0;
        grid_[left][j].edge_geometry.F_s = 0;
        grid_[left][j].edge_geometry.F_n = 0;
        grid_[left][j].edge_geometry.F_w = 0;

        grid_[right][j].edge_geometry.F_e = 0;
        grid_[right][j].edge_geometry.F_s = 0;
        grid_[right][j].edge_geometry.F_n = 0;
        grid_[right][j].edge_geometry.F_w = 0;

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

// Нужно проводить после каждой итерации алгоритма для соблюдения инвариантости.
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


void Grid::Cell::init_flows(Cell E, Cell W, Cell N, Cell S)
{
    edge_geometry.F_e = edge_geometry.S_z_e / 2 * (fields.u_z_P + E.fields.u_z_P) + edge_geometry.S_r_e / 2 * (fields.u_r_P + E.fields.u_r_P);
    edge_geometry.F_w = edge_geometry.S_z_w / 2 * (fields.u_z_P + W.fields.u_z_P) + edge_geometry.S_r_w / 2 * (fields.u_r_P + W.fields.u_r_P);
    edge_geometry.F_s = edge_geometry.S_z_s / 2 * (fields.u_z_P + S.fields.u_z_P) + edge_geometry.S_r_s / 2 * (fields.u_r_P + S.fields.u_r_P);
    edge_geometry.F_n = edge_geometry.S_z_n / 2 * (fields.u_z_P + N.fields.u_z_P) + edge_geometry.S_r_n / 2 * (fields.u_r_P + N.fields.u_r_P);
}

