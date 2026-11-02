#pragma once
#include <utility>

#include "../model_config.hpp"



class Grid
{
public:
    struct Cell
    {
        double theta_0, theta_1;
        double R_0, R_1;

        Cell(double theta_0, double theta_1, double R_0, double R_1);

        double R_P, theta_P;
        double r_P, d_P, V_P, sigma_P;
        double A_e, A_w, A_n, A_s;

        Vector2D S_e, S_w, S_n, S_s;

        bool is_boundary = false;
    };


    Grid(ModelConfig config);

    void init();

    std::pair<int, int> get_grid_dimensions();

private:

    void init_cells();
    void init_boundary_conditions(ModelConfig config);

    std::vector<std::vector<Cell>> grid_;

    int h_g, w_g;
}