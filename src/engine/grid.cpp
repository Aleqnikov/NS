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

Grid::Cell::Cell(double theta_0, double theta_1, double R_0, double R_1)
    : theta_0(theta_0), theta_1(theta_1), R_0(R_0), R_1(R_1)
{
    R_P = (R_0 - R_1) / 2;
    theta_P = (theta_1 - theta_0) / 2;
    double delta_theta = theta_1 - theta_0;
}
