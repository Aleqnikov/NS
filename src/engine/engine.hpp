#pragma once
#include "model_config.hpp"

struct Vector2D
{
    double x, y;
};

using PressureSpace = std::vector<Vector2D>;
using VelocitySpace = std::vector<Vector2D>;
using ProgramOutput = std::pair<PressureSpace, VelocitySpace>;




class Engine
{
public:
    Engine(ModelConfig config);
    void run();
    ProgramOutput get_results();


private:
    ModelConfig config_;
}