#include "engine.hpp"



Engine::Engine(ModelConfig config) : config_(config)
{
    grid_ = grid(config_);



}