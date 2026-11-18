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

#include "Model.hpp"
#include "Parameters.hpp"
#include "types.hpp"
#include "Grid.hpp"
#include "InitialCondition.hpp"
#include "ResultIO.hpp"
#include "RunFolder.hpp"
#include "enums.hpp"

// Prints one line about the progress after step_number frames of this run have been computed
// (counted from initial_time, so 1..frame_num, not the absolute frame number in the file name):
// the share of frames done, the estimated remaining time, the time the last frame took (elapsed, in ms)
// and the integration steps it took (accepted, and attempts rejected by the adaptivity).
void print_progress(Index step_number, const Parameters& parameters, int elapsed, long long accepted, long long rejected)
{
    double remaining_time = elapsed / 1000.0 * (parameters.frame_num - step_number);

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
              << std::setw(3) << elapsed_milisec << "ms";
    // After "Elapsed time", so the patterns of `cluster.sh status` still match the line.
    std::cout << "    Integration steps: " << accepted << " (rejected: " << rejected << ")" << std::endl;
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
    const Grid grid = parameters.grid();

    RunFolder run_folder(result_path);
    if (!run_folder.create())
        return EXIT_FAILURE;
    run_folder.save_setup(parameters, input_config_path);
    run_folder.write_run_start_info(run_start);

    Model model(parameters);
    Vector u(model.getDegreesOfFreedom());

    if (parameters.init_cond_from_file){
        // The path in the config is relative to the calculations/ folder of this run.
        const std::filesystem::path init_cond_path = run_folder.calculations() / parameters.init_cond_file_path;
        std::cout << "Setting initial condition from file: " << init_cond_path << std::endl;
        if (!set_init_cond_from_file(u, init_cond_path, grid, parameters.init_cond_scaling)){
            std::cerr << "Error: could not set the initial condition from the file " << init_cond_path
                      << " (init_cond_file_path = " << parameters.init_cond_file_path
                      << ", resolved relative to " << run_folder.calculations() << "). Exiting." << std::endl;
            return EXIT_FAILURE;
        }
    }
    else{
        std::cout << "Setting initial condition from code: " << ICType_to_string(parameters.init_condition) << std::endl;
        InitialCondition init_cond(parameters);
        set_init_cond_manually(u, init_cond, grid);
    }
    // The time stepping keeps the boundary values consistent only if they are consistent at the start.
    set_init_cond_boundary(u, grid, parameters.model_params());
    std::cout << "Initial condition set." << std::endl;

    // Wall time of the time stepping, including writing the frames.
    const auto computation_start = std::chrono::steady_clock::now();

    auto rhs = [ = ] __cuda_callable__( const TNL::Containers::StaticArray< 2, Index >& ind, const VectorView& u, VectorView& fu ) {
        model.set_rhs_at(u, fu, ind.x(), ind.y());
    };

    auto boundary_xdir = [ = ] __cuda_callable__( const Index ind, VectorView& fu ) {
        model.apply_boundary_condition_xdir(ind, fu);
    };

    auto boundary_ydir = [ = ] __cuda_callable__( const Index ind, VectorView& fu ) {
        model.apply_boundary_condition_ydir(ind, fu);
    };

    // Number of RHS evaluations. TNL calls time_stepping on the host once per stage of every attempted
    // step, including the attempts the adaptivity rejects, which it does not count itself.
    long long rhs_evaluations = 0;

    auto time_stepping = [ =, &rhs_evaluations ]( const Real& t, const Real& tau, const VectorView& u, VectorView& fu )
    {
        ++rhs_evaluations;

        // iterate over inner points only
        TNL::Containers::StaticArray< 2, Index > begin = {1, 1};
        TNL::Containers::StaticArray< 2, Index > end = {grid.sizeX - 1, grid.sizeY - 1};
        TNL::Algorithms::parallelFor< Device >(begin, end, rhs, u, fu );

        // iterate over boundary points (after the interior, Neumann copies the fu of the neighbour)
        TNL::Algorithms::parallelFor< Device >(0, grid.sizeX, boundary_xdir, fu );
        TNL::Algorithms::parallelFor< Device >(0, grid.sizeY, boundary_ydir, fu );
    };

    ODESolver solver;
    solver.setTau( parameters.integrationTimeStep );
    solver.setTime( parameters.initial_time );

    // Frame numbers are absolute (time / timeStep), so a run continued from a file with a later
    // initial_time carries on with the numbering of the previous run.
    const Index first_frame = static_cast<Index>(round(parameters.initial_time / parameters.timeStep));
    write_frame(run_folder.calculations(), first_frame, u, grid);

    long long accepted_steps = 0;  // integration steps accepted over the whole run

    // The loop counts frames instead of comparing the solver time with final_time, so it always ends
    // after frame_num frames, whatever the rounding of the time (a float/double mix once made it run forever).
    for (Index i = 1; i <= parameters.frame_num; ++i) {
        // Computed from i (in double), not by adding timeStep to the previous stop time, so the rounding
        // errors do not add up; the last frame stops exactly at final_time.
        const Real stop_time = (i == parameters.frame_num)
            ? parameters.final_time
            : static_cast<Real>(static_cast<double>(parameters.initial_time)
                                + (static_cast<double>(parameters.final_time) - static_cast<double>(parameters.initial_time))
                                  * i / parameters.frame_num);
        solver.setStopTime( stop_time );

        auto stepStart = std::chrono::high_resolution_clock::now();
        const long long rhs_evaluations_before = rhs_evaluations;
        solver.solve( u, time_stepping );
        write_frame(run_folder.calculations(), first_frame + i, u, grid);
        auto stepEnd = std::chrono::high_resolution_clock::now();

        // solve() resets the iteration counter, so getIterations() is the steps accepted in this frame.
        // Every attempt (accepted or rejected) evaluates the RHS once per stage of the method.
        const long long frame_accepted = solver.getIterations();
        const long long frame_attempts = (rhs_evaluations - rhs_evaluations_before) / Method::getStages();
        accepted_steps += frame_accepted;

        int elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(stepEnd - stepStart).count();
        print_progress(i, parameters, elapsed, frame_accepted, frame_attempts - frame_accepted);
    }

    const double computation_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - computation_start).count();

    if (rhs_evaluations % Method::getStages() != 0)
        std::cerr << "Warning: " << rhs_evaluations << " RHS evaluations is not a multiple of the "
                  << Method::getStages() << " stages of the method, so the number of rejected steps is wrong." << std::endl;
    const long long rejected_steps = rhs_evaluations / Method::getStages() - accepted_steps;
    std::cout << "Integration steps: " << accepted_steps << " (rejected: " << rejected_steps << ")" << std::endl;
    // Before the computation time: that line marks a finished run (cluster.sh), so everything else is written by then.
    run_folder.append_integration_steps(accepted_steps, rejected_steps);

    std::cout << "Computation time: " << format_computation_time(computation_seconds) << std::endl;
    run_folder.append_computation_time(computation_seconds);
}