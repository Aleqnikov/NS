#include "engine.hpp"
#include <fstream>
#include <cmath>
#include <iomanip>

Engine::Engine(ModelConfig config, std::string output_file) : config_(config), computation_grid_(config)
{
    std::ofstream fout(output_file);

    fout << "# vtk DataFile Version 3.0\n";
    fout << "sphere 2D semicircle\n";
    fout << "ASCII\n";
    fout << "DATASET STRUCTURED_GRID\n";

    std::size_t n = computation_grid_.grid_.size();
    std::size_t m = computation_grid_.grid_[0].size();

    std::size_t nx = m - 1;
    std::size_t ny = n - 1;

    fout << "DIMENSIONS " << (m - 1) << " " << (n - 1) << " 1\n";
    fout << "POINTS " << (m - 1) * (n - 1) << " double\n";

    fout << std::scientific << std::setprecision(8);

    for (std::size_t i = 0; i < n - 1; i++)
    {
        for (std::size_t j = 0; j < m - 1; j++)
        {
            double theta = (j == 0) ? computation_grid_.grid_[1][1].geometry.theta_s : computation_grid_.grid_[1][j].geometry.theta_n;
            double R     = (i == 0) ? computation_grid_.grid_[1][1].geometry.R_w : computation_grid_.grid_[i][1].geometry.R_e;

            fout << R * std::cos(theta) << " " << R * std::sin(theta) << " 0\n";
        }
    }

    fout << "\nCELL_DATA " << (nx - 1) * (ny - 1) << "\n";

    fout << "SCALARS Pressure double 1\n";
    fout << "LOOKUP_TABLE default\n";
    for (std::size_t i = 1; i < n - 1; i++)
    {
        for (std::size_t j = 1; j < m - 1; j++)
        {
            fout << computation_grid_.grid_[i][j].fields.p_P << "\n";
        }
    }

    fout << "\nVECTORS Velocity double\n";
    for (std::size_t i = 1; i < n - 1; i++)
    {
        for (std::size_t j = 1; j < m - 1; j++)
        {
            double uz = computation_grid_.grid_[i][j].fields.u_z_P;
            double ur = computation_grid_.grid_[i][j].fields.u_r_P;

            fout << uz << " " << ur << " 0\n";
        }
    }

}