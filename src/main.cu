#include <cstdlib>
#include <cassert>
#include <iostream>
#include <vector>
#include <filesystem>
#include <json.hpp>
#include <string>

#include "Problem.hpp"
#include "Parameters.hpp"
#include "types.hpp"
#include "InitialCondition.hpp"

int main(int argc, char** argv) {
    Parameters parameters = Parameters::load("config/config.json");
    std::filesystem::path result_path = "Results";
    if (argc == 2)
        result_path = argv[1];

    std::filesystem::path info_path   = result_path / "info";
    std::filesystem::path setup_path  = info_path / "parameters.txt";
    std::filesystem::path calc_path   = result_path / "calculations";
    std::filesystem::path config_path = info_path / "config.json";
    std::filesystem::path latex_path = info_path / "parameters.tex";
    
    try {
        std::filesystem::create_directories(result_path);
        std::filesystem::create_directories(info_path);
        std::filesystem::create_directories(calc_path);
    }
    catch (const std::filesystem::filesystem_error& e) {
        std::cerr << "Error creating folders: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    parameters.save_human_readable(setup_path);
    parameters.save_copy_of_config("config/config.json", config_path);
    parameters.save_for_latex(latex_path);

    Problem problem(parameters);
    Vector u(problem.getDegreesOfFreedom());

    if (!parameters.init_cond_from_file ||
        !problem.set_init_cond_from_file(u, calc_path / parameters.init_cond_file_path)) {
        std::cout << "Setting initial condition by code." << std::endl;
        InitialCondition init_cond(parameters.init_condition, parameters.domain, parameters.sizeX, parameters.sizeY, parameters.ksi);
        problem.set_init_cond_manually(u, init_cond);
    }

    ODESolver solver;
    solver.setTau( parameters.integrationTimeStep );
    solver.setTime( parameters.initial_time );
    
    Index step_number = static_cast<Index>(round(solver.getTime() / parameters.timeStep));
    problem.writeSolution(0.0, step_number, u, calc_path);
    
    while( solver.getTime() < parameters.final_time ) {
        solver.setStopTime( TNL::min( solver.getTime() + parameters.timeStep, parameters.final_time ) );

        auto rhs = [ = ] __cuda_callable__( const TNL::Containers::StaticArray< 2, Index >& ind, const VectorView& u, VectorView& fu ) mutable {
            problem.set_rhs_at(u, fu, ind.x(), ind.y());
        };

        auto boundary_xdir = [ = ] __cuda_callable__( const Index ind, VectorView& u, VectorView& fu ) mutable {
            problem.apply_boundary_condition_xdir(ind , u, fu);
        };

        auto boundary_ydir = [ = ] __cuda_callable__( const Index ind, VectorView& u, VectorView& fu ) mutable {
            problem.apply_boundary_condition_ydir(ind , u, fu);
        };

        auto time_stepping = [ = ]( const Real& t, const Real& tau, const VectorView& u, VectorView& fu )
        {
            // iterate over inner points only
            TNL::Containers::StaticArray< 2, Index > begin = {1, 1};
            TNL::Containers::StaticArray< 2, Index > end = {parameters.sizeX - 1, parameters.sizeY - 1};
            TNL::Algorithms::parallelFor< Device >(begin, end, rhs, u, fu );

            // iterate over boundary points
            TNL::Algorithms::parallelFor< Device >(0, parameters.sizeX, boundary_xdir, u, fu );
            TNL::Algorithms::parallelFor< Device >(0, parameters.sizeY, boundary_ydir, u, fu );
        };
        auto stepStart = std::chrono::high_resolution_clock::now();
        solver.solve( u, time_stepping );
        auto stepEnd = std::chrono::high_resolution_clock::now();
        step_number = static_cast<Index>(round(solver.getTime() / parameters.timeStep));
        problem.writeSolution(solver.getTime(), step_number, u, calc_path);

        double time_per_step = (std::chrono::duration<double>(stepEnd - stepStart).count());

        double remaining_time = time_per_step * (parameters.frame_num * parameters.final_time
                                                 / (parameters.final_time - parameters.initial_time)
                                                 - step_number);
      
        int hours = static_cast<int>(remaining_time) / 3600;
        int minutes = (static_cast<int>(remaining_time) % 3600) / 60;
        int seconds = remaining_time - (hours * 3600 + minutes * 60);
        
        std::cout << "Steps completed: " << step_number << " / " << parameters.frame_num << " => " << std::fixed
                  << std::setprecision(2) << ( double ) step_number / ( double ) parameters.frame_num * 100.0 << "% ";
        std::cout << "     Time remaining: " 
                  << hours << "h " << minutes << "m " << seconds << "s"
                  << std::endl;
    }
}