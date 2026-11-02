#pragma once


struct ModelConfig
{
    double speed;
    double radius;
    double density;     // Значение от 0 до 1, где 1 -- 100 процентов плотная.

    double rho = 1.225;         // кг/м^3
    double nu =  1.46e-5;       // м^2/c
    double nut =  4.38e-5;      //  м^2/c


    ModelConfig(double speed, double radius, double density) :
    speed{speed}, radius{radius}, density{density} {};
};