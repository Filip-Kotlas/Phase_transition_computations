#include <cstdlib>
#include <cassert>
#include <iostream>
#include <vector>
#include <filesystem>

#include <json.hpp>

#include "ode-solve.hpp"
#include "Merson.hpp"
#include "RungeKutta.hpp"
#include "ODEProblem.hpp"
#include "ODESolver.hpp"
#include "ACEProblem.hpp"
#include "Parameters.hpp"

namespace fs = std::filesystem;

int main(int argc, char** argv) {
    Parameters parameters = Parameters::load("config/config.json");

    fs::path result_path = "Results";
    if (argc == 2)
        result_path = argv[1];

    fs::path info_path   = result_path / "info";
    fs::path setup_path  = info_path / "parameters.txt";
    fs::path calc_path   = result_path / "calculations";
    fs::path config_path = info_path / "config.json";

    try {
        fs::create_directories(result_path);
        fs::create_directories(info_path);
        fs::create_directories(calc_path);
    }
    catch (const fs::filesystem_error& e) {
        std::cerr << "Error creating folders: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    parameters.save_human_readable(setup_path);
    parameters.save_copy_of_config("config/config.json", config_path);

    ACEProblem problem(parameters.sizeX,
                       parameters.sizeY,
                       parameters.domain,
                       parameters.alpha,
                       parameters.beta,
                       parameters.par_a,
                       parameters.ksi,
                       parameters.model,
                       calc_path);

    RungeKutta runge_kutta_integrator;
    Merson merson_integrator;

    double* u = new double[problem.getDegreesOfFreedom()];
    if (!parameters.init_cond_from_file ||
        !problem.set_init_cond_from_file(u, calc_path / parameters.init_cond_file_path)) {
        std::cout << "Setting initial condition by code." << std::endl;
        problem.set_init_cond_manually(u, parameters.init_condition);
    }
    problem.writeSolution(0.0, 0, u);

    if (parameters.type == "Runge-Kutta") {
        if (!solve(parameters.initial_time,
                   parameters.final_time,
                   parameters.timeStep,
                   parameters.integrationTimeStep,
                   &problem,
                   &runge_kutta_integrator,
                   u)) {
            delete[] u;
            return EXIT_FAILURE;
        }
    }
    else if (parameters.type == "Merson") {
        if (!solve(parameters.initial_time,
                   parameters.final_time,
                   parameters.timeStep,
                   parameters.integrationTimeStep,
                   &problem,
                   &merson_integrator,
                   u)) {
            delete[] u;
            return EXIT_FAILURE;
        }
    }

    delete[] u;
    return EXIT_SUCCESS;
}
