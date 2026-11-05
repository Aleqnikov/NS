#include "engine.hpp"

#include <iostream>
#include <ostream>

#include "solver.hpp"

Engine::Engine(ModelConfig config) : config_(config), computation_grid_(config)
{
}