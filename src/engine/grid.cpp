#include "grid.hpp"


Grid::Grid(ModelConfig config)
{
    double height = config.radius * 3;
    double width = config.radius  * 5;

    int N_R_max = 1000, N_R_min = 100;
    int N_theta_max = 1500, N_theta_min = 150;

    int N_R_current  = N_R_min + config.density * (N_R_max - N_R_min);
    int N_theta_current = N_theta_min + config.density * (N_theta_max - N_theta_min);

    grid_.resize(N_R_current + 2, N_theta_current + 2);


    this->init_cells();
    this->init_boundary_conditions(ModelConfig config);
}


std::pair<int, int> Grid::get_grid_dimensions()
{
    return {h_g, w_g};
}

Grid::Cell::Cell(double theta_0, double theta_1, double R_e, double R_w)
    : theta_0(theta_0), theta_1(theta_1), R_0(R_0), R_1(R_1)
{
    R_P = (R_0 - R_1) / 2;
    theta_P = (theta_1 - theta_0) / 2;
    double delta_theta = theta_1 - theta_0;

    double delta_c = std::cos(theta_1) - std::cos(theta_0);


    V_P = (R_e * R_e * R_e - R_w * R_w * R_w) / 3 * delta_c;
    sigma_P = (R_e * R_e - R_w * R_w) / 2 * delta_theta;

    A_e = R_e * R_e * delta_c;
    A_w = R_w * R_w * delta_c;
    A_n = std::sin(theta_1) * (R_e * R_e - R_w * R_w) / 2;
    A_s = std::sin(theta_2) * (R_e * R_e - R_w * R_w) / 2;

    S_n = {- std::sin(theta_n) * A_n, std::cos(theta_n) * A_n};
    S_s = {- std::sin(theta_s) * A_s, std::cos(theta_s) * A_s};

}
