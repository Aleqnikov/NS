#pragma once
#include <utility>

#include "../model_config.hpp"

struct SA_params
{
	double xi_P;
	double f_v1_P;
	double f_v2_P;
	double nu_t_P;
	double nu_eff_P;
	double gamma_P;
	double Omega_P;
	double S_tilda_P;
	double r_sa_P;
	double g_P;
	double f_w_P;
};

struct Simple
{
	double d_z_P;
	double d_r_P;
	double Theta_f;
	double p_stroke_P;
	double F_stroke_f;
	double u_stroke_z_P;
	double u_stroke_r_P;
};

struct Flows
{
	double F_f;
	double G_z_P;
	double G_r_P;
	double rho_z_f;
	double rho_r_f;
	double D_e, D_w, D_n, D_s;
	double a_P, a_E, a_W, a_N, a_S;

	double b_P;
}

class Grid
{
public:
    struct Cell
    {
        double theta_0, theta_1;
        double R_0, R_1;

        Cell(double theta_0, double theta_1, double R_0, double R_1);
        Cell(bool is_boundary);

        double R_P, theta_P;
        double r_P, d_P, V_P, sigma_P;
        double A_e, A_w, A_n, A_s;

        Vector2D S_e, S_w, S_n, S_s;

        bool is_boundary = false;

        // параметры считаемые динамически

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