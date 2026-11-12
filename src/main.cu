#include <cstdlib>
#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <ctime>
#include <vector>
#include <filesystem>
#include <json.hpp>
#include <string>

#include "Problem.hpp"
#include "Parameters.hpp"
#include "types.hpp"
#include "InitialCondition.hpp"
#include "ResultIO.hpp"
#include "RunFolder.hpp"
#include "enums.hpp"

// Prints one line about the progress after the frame step_number has been written:
// the share of frames done, the estimated remaining time and the time the last frame took (elapsed, in ms).
void print_progress(Index step_number, const Parameters& parameters, int elapsed)
{
    double remaining_time = elapsed / 1000.0 * (parameters.frame_num * parameters.final_time
                                             / (parameters.final_time - parameters.initial_time)
                                             - step_number);

    int remain_hours = static_cast<int>(remaining_time) / 3600;
    int remain_minutes = (static_cast<int>(remaining_time) % 3600) / 60;
    int remain_seconds = remaining_time - (remain_hours * 3600 + remain_minutes * 60);

    int elapsed_hours = elapsed / 3600000;
    int elapsed_minutes = (elapsed % 3600000) / 60000;
    int elapsed_seconds = (elapsed % 60000) / 1000;
    int elapsed_milisec = elapsed % 1000;

    std::cout << "Steps completed: " << std::setw(2) << step_number << " / " << std::setw(2) << parameters.frame_num << " => " << std::setw(5) << std::fixed
              << std::setprecision(2) << ( double ) step_number / ( double ) parameters.frame_num * 100.0 << "% ";
    std::cout << "    Time remaining: "
              << remain_hours << "h "
              << std::setw(2) << remain_minutes << "m "
              << std::setw(2) << remain_seconds << "s";
    std::cout << "    Elapsed time: "
              << elapsed_hours << "h "
              << std::setw(2) << elapsed_minutes << "m "
              << std::setw(2) << elapsed_seconds << "s "
              << std::setw(3) << elapsed_milisec << "ms" << std::endl;
}

int main(int argc, char** argv) {
    const std::time_t run_start = std::time(nullptr);

    // usage: Phase-field [result_folder [config_file]]
    std::filesystem::path result_path = "Results";
    if (argc >= 2)
        result_path = argv[1];
    std::filesystem::path input_config_path = "config/config.json";
    if (argc >= 3)
        input_config_path = argv[2];

    Parameters parameters = Parameters::load(input_config_path);

    RunFolder run_folder(result_path);
    if (!run_folder.create())
        return EXIT_FAILURE;
    run_folder.save_setup(parameters, input_config_path);
    run_folder.write_run_start_info(run_start);

    Problem problem(parameters);
    Vector u(problem.getDegreesOfFreedom());

    if (parameters.init_cond_from_file){
        // The path in the config is relative to the calculations/ folder of this run.
        const std::filesystem::path init_cond_path = run_folder.calculations() / parameters.init_cond_file_path;
        std::cout << "Setting initial condition from file: " << init_cond_path << std::endl;
        if (!set_init_cond_from_file(u, init_cond_path, parameters.sizeX, parameters.sizeY, parameters.init_cond_scaling)){
            std::cerr << "Error: could not set the initial condition from the file " << init_cond_path
                      << " (init_cond_file_path = " << parameters.init_cond_file_path
                      << ", resolved relative to " << run_folder.calculations() << "). Exiting." << std::endl;
            return EXIT_FAILURE;
        }
    }
    else{
        std::cout << "Setting initial condition from code: " << ICType_to_string(parameters.init_condition) << std::endl;
        InitialCondition init_cond(parameters);
        set_init_cond_manually(u, init_cond, parameters.sizeX, parameters.sizeY);
    }
    std::cout << "Initial condition set." << std::endl;

    // Wall time of the time stepping, including writing the frames.
    const auto computation_start = std::chrono::steady_clock::now();

    ODESolver solver;
    solver.setTau( parameters.integrationTimeStep );
    solver.setTime( parameters.initial_time );
    
    Index step_number = static_cast<Index>(round(solver.getTime() / parameters.timeStep));
    write_frame(run_folder.calculations(), step_number, u, parameters.sizeX, parameters.sizeY);

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
        step_number = static_cast<Index>(round(solver.getTime() / parameters.timeStep));
        write_frame(run_folder.calculations(), step_number, u, parameters.sizeX, parameters.sizeY);
        auto stepEnd = std::chrono::high_resolution_clock::now();

        int elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(stepEnd - stepStart).count();
        print_progress(step_number, parameters, elapsed);
    }

    const double computation_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - computation_start).count();
    std::cout << "Computation time: " << format_computation_time(computation_seconds) << std::endl;
    run_folder.append_computation_time(computation_seconds);
}