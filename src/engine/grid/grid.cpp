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
    int N_theta_max = 2000, N_theta_min = 100;

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
    this->init_flows();
    this->update_sa_parameters();

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
    double du_zd_theta = (N.fields.u_z_P - S.fields.u_z_P) /(N.geometry.theta_P - S.geometry.theta_P);


    double du_zdr = std::sin(geometry.theta_P) * du_zd_R +
        std::cos(geometry.theta_P)/geometry.R_P * du_zd_theta;


    double du_rd_R =  (E.fields.u_r_P - W.fields.u_r_P)/(E.geometry.R_P - W.geometry.R_P);
    double du_rd_theta = (N.fields.u_r_P - S.fields.u_r_P) /(N.geometry.theta_P - S.geometry.theta_P);

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

        grid_[i][border].geometry.R_P = grid_[i][border + 1].geometry.R_P;
        grid_[i][border].geometry.theta_P = 2.0 * grid_[i][border + 1].geometry.theta_s  -grid_[i][border + 1].geometry.theta_P;
        grid_[i][border].geometry.r_P = grid_[i][border].geometry.R_P * std::sin(grid_[i][border].geometry.theta_P);
        grid_[i][upper] = Cell(true);
        grid_[i][upper].fields.nut_P = grid_[i][upper - 1].fields.nut_P;
        grid_[i][upper].fields.p_P = grid_[i][upper - 1].fields.p_P;
        grid_[i][upper].fields.u_r_P = -grid_[i][upper - 1].fields.u_r_P;
        grid_[i][upper].fields.u_z_P = grid_[i][upper - 1].fields.u_z_P;
        grid_[i][upper].geometry.R_P = grid_[i][upper - 1].geometry.R_P;
        grid_[i][upper].geometry.theta_P = 2.0 * grid_[i][upper - 1].geometry.theta_n - grid_[i][upper - 1].geometry.theta_P;
        grid_[i][upper].geometry.r_P = grid_[i][upper].geometry.R_P * std::sin(grid_[i][upper].geometry.theta_P);
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

        grid_[left][j].geometry.R_P = 2.0 * grid_[left + 1][j].geometry.R_w - grid_[left + 1][j].geometry.R_P;
        grid_[left][j].geometry.theta_P = grid_[left + 1][j].geometry.theta_P;
        grid_[left][j].geometry.r_P = grid_[left][j].geometry.R_P * std::sin(grid_[left][j].geometry.theta_P);
        grid_[right][j].geometry.R_P = 2.0 * grid_[right - 1][j].geometry.R_e - grid_[right - 1][j].geometry.R_P;
        grid_[right][j].geometry.theta_P = grid_[right - 1][j].geometry.theta_P;
        grid_[right][j].geometry.r_P = grid_[right][j].geometry.R_P * std::sin(grid_[right][j].geometry.theta_P);


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

void Grid::init_flows()
{
    for (std::size_t i = 1; i < grid_.size() - 1; ++i)
    {
        for (std::size_t j = 1; j < grid_[i].size() - 1; ++j)
        {
            grid_[i][j].init_flows(grid_[i + 1][j], grid_[i - 1][j], grid_[i][j + 1], grid_[i][j - 1]);
        }
    }
}

void Grid::update_sa_parameters()
{
    for (std::size_t i = 1; i < grid_.size() - 1; ++i)
    {
        for (std::size_t j = 1; j < grid_[i].size() - 1; ++j)
        {
            grid_[i][j].update_SA_parameters(grid_[i + 1][j], grid_[i - 1][j], grid_[i][j + 1], grid_[i][j - 1], config_);
        }
    }
}

std::size_t pair_idx_to_uno(std::size_t i, std::size_t j, std::size_t n)
{
    // Нужно учесть то, что у нас как бы ни i ни j не должны быть граничными и тп.
    return i * n + j;
}

void Grid::step_speeds_update()
{

    std::vector<std::vector<double>> A(grid_.size() * grid_.size(), grid_.size() * grid_.size());

    std::vector<double> b(grid_.size() * grid_.size());

    for (std::size_t i = 1; i < grid_.size() - 1; ++i)
    {
        for (std::size_t j = 1; j < grid_[i].size() - 1; ++j)
        {

            Cell E = grid_[i + 1][j];
            Cell W = grid_[i - 1][j];
            Cell N = grid_[i][j + 1];
            Cell S = grid_[i][j - 1];

            Cell P = grid_[i][j];

            double alpha_z_e = 2 * P.edge_geometry.S_z_e * std::cos(P.geometry.theta_P) + P.edge_geometry.S_r_e * std::sin(P.geometry.theta_P);
            double beta_z_e = - 2 * P.edge_geometry.S_z_e * std::sin(P.geometry.theta_P) + P.edge_geometry.S_r_e * std::cos(P.geometry.theta_P);

            double alpha_z_w = 2 * P.edge_geometry.S_z_w * std::cos(P.geometry.theta_P) + P.edge_geometry.S_r_w * std::sin(P.geometry.theta_P);
            double beta_z_w = - 2 * P.edge_geometry.S_z_w * std::sin(P.geometry.theta_P) + P.edge_geometry.S_r_w * std::cos(P.geometry.theta_P);

            double alpha_z_n = 2 * P.edge_geometry.S_z_n * std::cos(P.geometry.theta_n) + P.edge_geometry.S_r_n * std::sin(P.geometry.theta_n);
            double beta_z_n = - 2 * P.edge_geometry.S_z_n * std::sin(P.geometry.theta_n) + P.edge_geometry.S_r_n * std::cos(P.geometry.theta_n);

            double alpha_z_s = 2 * P.edge_geometry.S_z_s * std::cos(P.geometry.theta_s) + P.edge_geometry.S_r_s * std::sin(P.geometry.theta_s);
            double beta_z_s = - 2 * P.edge_geometry.S_z_s * std::sin(P.geometry.theta_s) + P.edge_geometry.S_r_s * std::cos(P.geometry.theta_s);

            double nu_e = (E.is_boundary ? P.sa_params.nu_eff_P : (E.sa_params.nu_eff_P + P.sa_params.nu_eff_P) / 2);
            double nu_w = (W.is_boundary ? P.sa_params.nu_eff_P : (W.sa_params.nu_eff_P + P.sa_params.nu_eff_P) / 2);
            double nu_s = (S.is_boundary ? P.sa_params.nu_eff_P : (S.sa_params.nu_eff_P + P.sa_params.nu_eff_P) / 2);
            double nu_n = (N.is_boundary ? P.sa_params.nu_eff_P : (N.sa_params.nu_eff_P + P.sa_params.nu_eff_P) / 2);

            double D_e = (nu_e*alpha_z_e)/(E.geometry.R_P - P.geometry.R_P);
            double D_w = (nu_w * alpha_z_w)/(P.geometry.R_P - W.geometry.R_P);
            double D_n = (nu_n * beta_z_n)/(P.geometry.R_P * (N.geometry.theta_P - P.geometry.theta_P));
            double D_s = (nu_s * beta_z_s)/(P.geometry.R_P * (P.geometry.theta_P - S.geometry.theta_P));

            double a_E = D_e + std::max(-P.edge_geometry.F_e, 0.);
            double a_W = D_w + std::max(P.edge_geometry.F_w, 0.);
            double a_N = D_n + std::max(-P.edge_geometry.F_n, 0.);
            double a_S = D_s + std::max(P.edge_geometry.F_s, 0.);
            double a_P = a_E + a_W + a_N + a_S;


            double G_z_P = ((E.fields.p_P + P.fields.p_P) / 2 - P.fields.p_P) * P.edge_geometry.S_z_e -
                           ((W.fields.p_P + P.fields.p_P) / 2 - P.fields.p_P) * P.edge_geometry.S_z_w +
                           ((N.fields.p_P + P.fields.p_P) / 2 - P.fields.p_P) * P.edge_geometry.S_z_n -
                           ((S.fields.p_P + P.fields.p_P) / 2 - P.fields.p_P) * P.edge_geometry.S_z_s;


            // считаем ро, боже, какой ужас. почему так много формул..... это же вообще не поддерживаемо...

            double du_zdR_P = (E.fields.u_z_P - W.fields.u_z_P) / (E.geometry.R_P - W.geometry.R_P);
            double du_zdtheta_P = (N.fields.u_z_P - S.fields.u_z_P)/(N.geometry.theta_P - S.geometry.theta_P);
            double du_rdR_p = (E.fields.u_r_P - W.fields.u_r_P) / (E.geometry.R_P - W.geometry.R_P);
            double du_rdtheta_P = (N.fields.u_r_P - S.fields.u_r_P)/(N.geometry.theta_P - S.geometry.theta_P);
            double du_rdz_P = std::cos(P.geometry.theta_P) * du_rdR_p - std::sin(P.geometry.theta_P) / P.geometry.R_P * du_rdtheta_P;

            double du_zdR_E = E.is_boundary ? du_zdR_P : (grid_[i+2][j].fields.u_z_P - grid_[i][j].fields.u_z_P) / (grid_[i+2][j].geometry.R_P - grid_[i][j].geometry.R_P);
            double du_zdR_W = W.is_boundary ? du_zdR_P : (grid_[i][j].fields.u_z_P - grid_[i-2][j].fields.u_z_P) / (grid_[i][j].geometry.R_P - grid_[i-2][j].geometry.R_P);
            double du_zdR_N = N.is_boundary ? du_zdR_P : (grid_[i+1][j+1].fields.u_z_P - grid_[i-1][j+1].fields.u_z_P) / (grid_[i+1][j+1].geometry.R_P - grid_[i-1][j+1].geometry.R_P);
            double du_zdR_S = S.is_boundary ? du_zdR_P : (grid_[i+1][j-1].fields.u_z_P - grid_[i-1][j-1].fields.u_z_P) / (grid_[i+1][j-1].geometry.R_P - grid_[i-1][j-1].geometry.R_P);

            double du_zdtheta_E = E.is_boundary ? du_zdtheta_P : (grid_[i+1][j+1].fields.u_z_P - grid_[i+1][j-1].fields.u_z_P) / (grid_[i+1][j+1].geometry.theta_P - grid_[i+1][j-1].geometry.theta_P);
            double du_zdtheta_W = W.is_boundary ? du_zdtheta_P : (grid_[i-1][j+1].fields.u_z_P - grid_[i-1][j-1].fields.u_z_P) / (grid_[i-1][j+1].geometry.theta_P - grid_[i-1][j-1].geometry.theta_P);
            double du_zdtheta_N = N.is_boundary ? du_zdtheta_P : (grid_[i][j+2].fields.u_z_P - grid_[i][j].fields.u_z_P)/ (grid_[i][j+2].geometry.theta_P - grid_[i][j].geometry.theta_P);
            double du_zdtheta_S = S.is_boundary ? du_zdtheta_P : (grid_[i][j].fields.u_z_P - grid_[i][j-2].fields.u_z_P) / (grid_[i][j].geometry.theta_P - grid_[i][j-2].geometry.theta_P);


            double du_zdR_e = 0.5 * (du_zdR_P + du_zdR_E);
            double du_zdR_w = 0.5 * (du_zdR_P + du_zdR_W);
            double du_zdR_n = 0.5 * (du_zdR_P + du_zdR_N);
            double du_zdR_s = 0.5 * (du_zdR_P + du_zdR_S);


            double du_rdR_E = E.is_boundary ? du_rdR_p : (grid_[i+2][j].fields.u_r_P - grid_[i][j].fields.u_r_P) / (grid_[i+2][j].geometry.R_P - grid_[i][j].geometry.R_P);
            double du_rdR_W = W.is_boundary ? du_rdR_p : (grid_[i][j].fields.u_r_P - grid_[i-2][j].fields.u_r_P) / (grid_[i][j].geometry.R_P - grid_[i-2][j].geometry.R_P);
            double du_rdR_N = N.is_boundary ? du_rdR_p : (grid_[i+1][j+1].fields.u_r_P - grid_[i-1][j+1].fields.u_r_P) / (grid_[i+1][j+1].geometry.R_P - grid_[i-1][j+1].geometry.R_P);
            double du_rdR_S = S.is_boundary ? du_rdR_p : (grid_[i+1][j-1].fields.u_r_P - grid_[i-1][j-1].fields.u_r_P) / (grid_[i+1][j-1].geometry.R_P - grid_[i-1][j-1].geometry.R_P);

            double du_rdtheta_E = E.is_boundary ? du_rdtheta_P : (grid_[i+1][j+1].fields.u_r_P - grid_[i+1][j-1].fields.u_r_P) / (grid_[i+1][j+1].geometry.theta_P - grid_[i+1][j-1].geometry.theta_P);
            double du_rdtheta_W = W.is_boundary ? du_rdtheta_P : (grid_[i-1][j+1].fields.u_r_P - grid_[i-1][j-1].fields.u_r_P) / (grid_[i-1][j+1].geometry.theta_P - grid_[i-1][j-1].geometry.theta_P);
            double du_rdtheta_N = N.is_boundary ? du_rdtheta_P : (grid_[i][j+2].fields.u_r_P - grid_[i][j].fields.u_r_P) / (grid_[i][j+2].geometry.theta_P - grid_[i][j].geometry.theta_P);
            double du_rdtheta_S = S.is_boundary ? du_rdtheta_P : (grid_[i][j].fields.u_r_P - grid_[i][j-2].fields.u_r_P) / (grid_[i][j].geometry.theta_P - grid_[i][j-2].geometry.theta_P);

            double du_rdz_E = E.is_boundary ? du_rdz_P : std::cos(E.geometry.theta_P)*du_rdR_E - std::sin(E.geometry.theta_P)/E.geometry.R_P*du_rdtheta_E;
            double du_rdz_W = W.is_boundary ? du_rdz_P : std::cos(W.geometry.theta_P)*du_rdR_W - std::sin(W.geometry.theta_P)/W.geometry.R_P*du_rdtheta_W;
            double du_rdz_N = N.is_boundary ? du_rdz_P : std::cos(N.geometry.theta_P)*du_rdR_N - std::sin(N.geometry.theta_P)/N.geometry.R_P*du_rdtheta_N;
            double du_rdz_S = S.is_boundary ? du_rdz_P : std::cos(S.geometry.theta_P)*du_rdR_S - std::sin(S.geometry.theta_P)/S.geometry.R_P*du_rdtheta_S;


            double du_zdtheta_e = E.is_boundary ? du_zdtheta_P : 0.5*(du_zdtheta_P + du_zdtheta_E);
            double du_zdtheta_w = W.is_boundary ? du_zdtheta_P : 0.5*(du_zdtheta_P + du_zdtheta_W);

            double du_rdz_e = E.is_boundary ? du_rdz_P : 0.5*(du_rdz_P + du_rdz_E);
            double du_rdz_w = W.is_boundary ? du_rdz_P : 0.5*(du_rdz_P + du_rdz_W);
            double du_rdz_n = N.is_boundary ? du_rdz_P : 0.5*(du_rdz_P + du_rdz_N);
            double du_rdz_s = S.is_boundary ? du_rdz_P : 0.5*(du_rdz_P + du_rdz_S);


            double rho_z_e = nu_e * (beta_z_e * 1/E.geometry.R_P * (du_zdtheta_e) + P.edge_geometry.S_r_e * du_rdz_e);
            double rho_z_w = nu_w * (beta_z_e * 1/W.geometry.R_P * (du_zdtheta_w) + P.edge_geometry.S_r_w * du_rdz_w);

            double rho_z_n = nu_n * (alpha_z_n * du_zdR_n + P.edge_geometry.S_r_n * du_rdz_n);
            double rho_z_s = nu_s * (alpha_z_s * du_zdR_s + P.edge_geometry.S_r_s * du_rdz_s);

            double b_P = -G_z_P + rho_z_e - rho_z_w + rho_z_n - rho_z_s;

            std::size_t new_i = pair_idx_to_uno(i, j, grid_.size());

            b[new_i] = b_P;

            A[new_i][new_i] = a_P;

            std::size_t new_e = pair_idx_to_uno(i + 1, j, grid_.size());
            std::size_t new_w = pair_idx_to_uno(i - 1, j, grid_.size());
            std::size_t new_n = pair_idx_to_uno(i, j + 1, grid_.size());
            std::size_t new_s = pair_idx_to_uno(i, j - 1, grid_.size());

            if (!E.is_boundary) {A[new_i][new_e] = a_E;}
            if (!W.is_boundary) {A[new_i][new_w] = a_W;}
            if (!N.is_boundary) {A[new_i][new_n] = a_N;}
            if (!S.is_boundary) {A[new_i][new_s] = a_S;}

            grid_[i][j].update_SA_parameters(grid_[i + 1][j], grid_[i - 1][j], grid_[i][j + 1], grid_[i][j - 1], config_);
        }
    }
}


