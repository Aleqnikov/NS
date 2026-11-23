#pragma once


struct ModelConfig
{
    double speed;
    double radius;
    double density;     // Значение от 0 до 1, где 1 -- 100 процентов плотная.

    double rho = 1.225;         // кг/м^3
    double nu =  1.46e-5;       // м^2/c
    double nut =  3 * nu;      //  м^2/c


    // константы для SA
    double c_v_1 = 7.1;
    double c_b_1 = 0.1355;
    double c_b_2 = 0.622;
    double k = 0.41;
    double c_w_2 = 0.3;
    double c_w_3 = 2;
    double sigma  = 2. / 3.;
	double c_w_1 = c_b_1 / (k * k) + (1 + c_b_2)/(sigma);


    // Параметры релаксации для подсчета давления
    double alpha_p = 0.3;
    double alpha_u = 0.7;

	double alpha_sa = 0.7;


    ModelConfig(double speed, double radius, double density) :
    speed{speed}, radius{radius}, density{density} {};
};