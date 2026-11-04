#pragma once


struct ModelConfig
{
    double speed;
    double radius;
    double density;     // Значение от 0 до 1, где 1 -- 100 процентов плотная.

    double rho = 1.225;         // кг/м^3
    double nu =  1.46e-5;       // м^2/c
    double nut =  4.38e-5;      //  м^2/c


    // константы для SA
    double c_v_1 = 7.1;
    double c_b_1 = 0.1355;
    double c_b_2 = 0.622;
    double k = 0.41;
    double c_w_2 = 0.3;
    double c_w_3 = 2;
    double sigma  = 0.6666666;


    ModelConfig(double speed, double radius, double density) :
    speed{speed}, radius{radius}, density{density} {};
};