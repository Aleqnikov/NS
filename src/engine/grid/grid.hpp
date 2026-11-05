#pragma once
#include <utility>



#include "../../model_config.hpp"
#include <vector>

#include "grid_structures.hpp"




class Grid
{
public:
    struct Cell
    {
        double theta_0, theta_1;
        double R_0, R_1;

        Cell() = default;
        Cell(double theta_0, double theta_1, double R_0, double R_1, double R_s, ModelConfig config);
        Cell(double theta_0_, double theta_1_, double R_0_, double R_1_, double R_s, ModelConfig config, bool is_boundary_);
        Cell(bool is_boundary_);

        void update_SA_parameters(Cell E, Cell W, Cell N, Cell S, ModelConfig config);
        void init_flows(Cell E, Cell W, Cell N, Cell S);

        bool is_boundary = false;

       	Fields fields;
       	CellGeometry geometry;
		EdgeGeometry edge_geometry;
        SA_params sa_params;
    };

    Grid(ModelConfig config);

    // Проходит по всем ячейкам и обновляет параметры SA считаемые по формулам.
    void update_sa_parameters();

    void init_flows();
    void step_speeds_update();
    void print();

    std::pair<int, int> get_grid_dimensions();

    std::vector<std::vector<Cell>> grid_;

private:

    int h_g, w_g;

    ModelConfig config_;


    void init_boundary_cond();
    void update_boundary_cells();
};