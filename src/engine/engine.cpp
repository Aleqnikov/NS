#include "engine.hpp"

#include <iostream>
#include <ostream>

#include "grid/solver.hpp"

Engine::Engine(ModelConfig config) : config_(config), computation_grid_(config)
{
}