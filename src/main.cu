#include <cstdlib>
#include <cassert>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <chrono>
#include <ctime>
#include <vector>
#include <filesystem>
#include <json.hpp>
#include <TNL/Backend/Functions.h>
#include <TNL/Backend/DeviceInfo.h>
#include <string>
#include <unistd.h>

#include "Problem.hpp"
#include "Parameters.hpp"
#include "types.hpp"
#include "InitialCondition.hpp"
#include "enums.hpp"

// Name and compute capability of the GPU this process runs on, e.g. "NVIDIA GeForce RTX 3060 (compute capability 8.6)".
std::string gpu_description()
{
    try {
        const int device = TNL::Backend::getDevice();
        return TNL::Backend::getDeviceName(device) + " (compute capability "
               + std::to_string(TNL::Backend::getArchitectureMajor(device)) + "."
               + std::to_string(TNL::Backend::getArchitectureMinor(device)) + ")";
    }
    catch (const std::exception&) {
        // TNL throws when the CUDA call fails; the run itself does not depend on this information.
        return "unknown";
    }
}

// Short name of the computer this process runs on, e.g. "gp11" on the cluster
// (the system returns "gp11.fjfi.cvut.cz"; the domain is dropped).
std::string host_name()
{
    char name[256] = {};
    // The last byte stays '\0' even if the name is truncated.
    if (gethostname(name, sizeof(name) - 1) != 0)
        return "unknown";
    const std::string full(name);
    return full.substr(0, full.find('.'));
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

    // info/setup holds what the run was started with; the Python scripts add info/plots and info/boundary.
    std::filesystem::path setup_dir   = result_path / "info" / "setup";
    std::filesystem::path setup_path  = setup_dir / "parameters.txt";
    std::filesystem::path calc_path   = result_path / "calculations";
    std::filesystem::path config_path = setup_dir / "config.json";
    std::filesystem::path latex_path  = setup_dir / "parameters.tex";

    try {
        std::filesystem::create_directories(setup_dir);
        std::filesystem::create_directories(calc_path);
    }
    catch (const std::filesystem::filesystem_error& e) {
        std::cerr << "Error creating folders: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    parameters.save_human_readable(setup_path);
    parameters.save_copy_of_config(input_config_path, config_path);
    parameters.save_for_latex(latex_path);

    // Written at the start, so that a run that crashes or is killed still records them.
    // The computation time is appended at the end and "make run" appends the git commit after it.
    {
        std::ofstream setup_file(setup_path, std::ios::app);
        setup_file << "\nRun date: " << std::put_time(std::localtime(&run_start), "%d.%m.%Y %H:%M:%S") << '\n'
                   << "Computer: " << host_name() << '\n'
                   << "GPU: " << gpu_description() << '\n';
        if (!setup_file)
            std::cerr << "Warning: unable to write the run date, computer and GPU to " << setup_path << std::endl;
    }

    Problem problem(parameters);
    Vector u(problem.getDegreesOfFreedom());

    if (parameters.init_cond_from_file){
        // The path in the config is relative to the calculations/ folder of this run.
        const std::filesystem::path init_cond_path = calc_path / parameters.init_cond_file_path;
        std::cout << "Setting initial condition from file: " << init_cond_path << std::endl;
        if (!problem.set_init_cond_from_file(u, init_cond_path, parameters.init_cond_scaling)){
            std::cerr << "Error: could not set the initial condition from the file " << init_cond_path
                      << " (init_cond_file_path = " << parameters.init_cond_file_path
                      << ", resolved relative to " << calc_path << "). Exiting." << std::endl;
            return EXIT_FAILURE;
        }
    }
    else{
        std::cout << "Setting initial condition from code: " << ICType_to_string(parameters.init_condition) << std::endl;
        InitialCondition init_cond(parameters);
        problem.set_init_cond_manually(u, init_cond);
    }
    std::cout << "Initial condition set." << std::endl;

    // Wall time of the time stepping, including writing the frames.
    const auto computation_start = std::chrono::steady_clock::now();

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
        step_number = static_cast<Index>(round(solver.getTime() / parameters.timeStep));
        problem.writeSolution(solver.getTime(), step_number, u, calc_path);
        auto stepEnd = std::chrono::high_resolution_clock::now();

        int elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(stepEnd - stepStart).count();

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

    const double computation_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - computation_start).count();
    const long total_seconds = std::lround(computation_seconds);
    std::ostringstream computation_time;
    computation_time << total_seconds / 3600 << "h "
                     << std::setfill('0') << std::setw(2) << (total_seconds % 3600) / 60 << "m "
                     << std::setw(2) << total_seconds % 60 << "s"
                     << " (" << std::fixed << std::setprecision(1) << computation_seconds << " s)";

    std::cout << "Computation time: " << computation_time.str() << std::endl;

    std::ofstream setup_file(setup_path, std::ios::app);
    setup_file << "Computation time: " << computation_time.str() << '\n';
    if (!setup_file)
        std::cerr << "Warning: unable to write the computation time to " << setup_path << std::endl;
}