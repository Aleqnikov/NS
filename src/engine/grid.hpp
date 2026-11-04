#pragma once
#include <utility>

#include "../model_config.hpp"
#include <vector>


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

struct Fields
{
	double u_z_P, u_r_P;
	double p_P;
	double nut_P;
};

struct CellGeometry
{
    CellGeometry() = default;
	CellGeometry(double R_w, double R_e, double theta_s, double theta_n, double R_s);

	double R_w, R_e;
	double theta_s, theta_n;

	double R_P;
	double theta_P;
	double r_P;
	double d_P;
	double delta_theta;
	double delta_c;
	double V_P;
	double sigma_P;
};

struct EdgeGeometry
{
    EdgeGeometry() = default;
	EdgeGeometry(CellGeometry cell_data);

	double A_e, A_w, A_n, A_s;
	double S_z_e, S_r_e;
	double S_z_w, S_r_w;
	double S_z_n, S_r_n;
	double S_z_s, S_r_s;

	// Потоки на гранях
	double F_e, F_w, F_n, F_s;
};


class Grid
{
public:
    struct Cell
    {
        double theta_0, theta_1;
        double R_0, R_1;

        Cell() = default;
        Cell(double theta_0, double theta_1, double R_0, double R_1, double R_s);
        Cell(bool is_boundary);



        bool is_boundary = false;

       	Fields fields;
       	CellGeometry geometry;
		EdgeGeometry edge_geometry;
        SA_params sa_params;
    };

    Grid(ModelConfig config);

    void init();

    void print();


    std::pair<int, int> get_grid_dimensions();

private:


    std::vector<std::vector<Cell>> grid_;

    int h_g, w_g;
};