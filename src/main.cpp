#include <iostream>
#include "CLI11.hpp"


#include "model_config.hpp"
#include "engine.hpp"



int main(int argc, char *argv[])
{
    double speed;
    double radius;
    double mesh_density;

    std::string output_file;

    CLI::App app("NS simulation");

    app.add_option("-s,--speed", speed, "Speed of the incoming flow");
    app.add_option("-r,--radius", radius, "The radius of the sphere");
    app.add_option("-d,--density", mesh_density, "The density of the mesh around the sphere");
    app.add_option("-o,--output", mesh_density, "Output file");

    app.parse(argc, argv);

    ModelConfig model_config(speed, radius, mesh_density);

    Engine engine(model_config);
    engine.run();
    auto results = engine.get_results();

    VTKConvertor vtk_convertor(results);

    vtk_convertor.convert(output_file);


    return 1;
}