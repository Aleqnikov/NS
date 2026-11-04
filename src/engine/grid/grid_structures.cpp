#include "grid_structures.hpp"

#include <cmath>

CellGeometry::CellGeometry(double R_w, double R_e, double theta_s, double theta_n, double R_s)
{
    this->R_w = R_w;
    this->R_e = R_e;
    this->theta_s = theta_s;
    this->theta_n = theta_n;

    R_P = (R_w + R_e) / 2;
    theta_P = (theta_s  + theta_n) / 2;
    r_P  = R_P * std::sin(theta_P);
    d_P = R_P - R_s;
    delta_theta = theta_n - theta_s;
    delta_c = std::cos(theta_s) - std::cos(theta_n);
    V_P = (R_e * R_e * R_e - R_w * R_w * R_w) / 3  * delta_c;
    sigma_P  = (R_e * R_e - R_w * R_w ) / 2  * delta_theta;
}

EdgeGeometry::EdgeGeometry(CellGeometry gm)
{


    A_e = gm.R_e * gm.R_e * gm.delta_c;
    A_w = gm.R_w * gm.R_w *  gm.delta_c;
    A_n = std::sin(gm.theta_n) * (gm.R_e * gm.R_e - gm.R_w * gm.R_w )/ 2;
    A_s = std::sin(gm.theta_s) * (gm.R_e * gm.R_e - gm.R_w * gm.R_w )/ 2;

    S_z_e = gm.R_e * gm.R_e / 2 * (std::sin(gm.theta_n) * std::sin(gm.theta_n) - std::sin(gm.theta_s) * std::sin(gm.theta_s));
    S_r_e = gm.R_e * gm.R_e / 2  * (gm.delta_theta - (std::sin(2 * gm.theta_n)  - std::sin(2 * gm.theta_s)) / 2);

    S_z_w = gm.R_w * gm.R_w / 2 * (std::sin(gm.theta_n) * std::sin(gm.theta_n) - std::sin(gm.theta_s) * std::sin(gm.theta_s));
    S_r_w = gm.R_w * gm.R_w / 2  * (gm.delta_theta - (std::sin(2 * gm.theta_n)  - std::sin(2 * gm.theta_s)) / 2);

    S_z_n = - A_n * std::sin(gm.theta_n);
    S_r_n = A_n * std::cos(gm.theta_n);

    S_z_s = - A_s * std::sin(gm.theta_s);
    S_r_s = A_s *  std::cos(gm.theta_s);

}


Fields::Fields(double speed, double nu)
{
    u_z_P = speed;
    u_r_P = 0;
    p_P = 0;
    nut_P = 3 * nu;
}