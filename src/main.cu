#include <cstdlib>
#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <ctime>
#include <vector>
#include <filesystem>
#include "external/json.hpp"
#include <limits>
#include <sstream>
#include <string>
#include <stdexcept>
#include <type_traits>

#include <TNL/Algorithms/reduce.h>

#include "models/AllenCahn.hpp"
#include "models/ZrNbChoud2013.hpp"
#include "models/ZrHMa2002.hpp"
#include "config/Parameters.hpp"
#include "core/types.hpp"
#include "core/Grid.hpp"
#include "core/warnings.hpp"
#include "initial_condition/InitialCondition.hpp"
#include "io/ResultIO.hpp"
#include "io/RunFolder.hpp"
#include "config/enums.hpp"

// Exit code of a run that started computing but stopped because the solution broke down (NaN or inf
// in a frame, or the solver stuck, see SolverStuck). The frames up to then are written (with NaN/inf in
// a frame, including the broken one), so `make run` still plots them
// (the Makefile checks for this value); any other failure, including a frame that cannot be
// written, exits with EXIT_FAILURE.
constexpr int EXIT_RUN_FAILED = 2;

// The solver counts as stuck after this many attempts in a row at the same time. With a NaN in the
// error estimate, Kutta-Merson neither accepts the step nor shrinks it (all comparisons with NaN are
// false) and retries it forever; with a collapsed step the time stops advancing. A healthy run rejects
// a step only rarely (1-3 times in thousands of steps), and each rejection shrinks the step at least
// by 0.8, so 100 rejections in a row would mean a step 0.8^100 ~ 2e-10 times smaller.
constexpr int max_attempts_at_same_time = 100;

// Thrown from time_stepping when the solver is stuck (it is the only way out of TNL's loop); caught in run.
struct SolverStuck : std::runtime_error {
    using std::runtime_error::runtime_error;
};

// Reports why the computation stops early: prints the reason and records it in parameters.txt
// ("Run failed: <reason>"), so the run does not look finished (no "Computation time" follows).
void record_failure(const RunFolder& run_folder, const std::string& reason)
{
    std::cerr << "Error: " << reason << ". Stopping the computation." << std::endl;
    run_folder.append_failure(reason);
}

// Number of values of u with indices in [begin, end) that are NaN or infinite (counted on the GPU).
Index count_non_finite(const Vector& u, Index begin, Index end)
{
    const auto view = u.getConstView();
    auto non_finite = [ = ] __cuda_callable__( Index k ) -> Index {
        // std::isfinite works in device code with nvcc (checked on the cluster) and, unlike the bare
        // CUDA isfinite, is also known to the host-only IntelliSense parser.
        return std::isfinite(view[k]) ? 0 : 1;
    };
    return TNL::Algorithms::reduce< Device >(begin, end, non_finite, TNL::Plus{});
}

// Checks that frame number `frame` (just written, at time t) contains no NaN or inf. If it does,
// prints the reason, records it in parameters.txt ("Run failed: ...") and returns false.
bool check_frame_finite(const Vector& u, const Grid& grid, Index frame, Real t, const RunFolder& run_folder)
{
    const Index N = grid.num_points();
    const Index bad_phase = count_non_finite(u, 0, N);
    const Index bad_conc = count_non_finite(u, N, 2 * N);
    if (bad_phase == 0 && bad_conc == 0)
        return true;

    std::ostringstream reason;
    reason << "NaN or inf in frame " << frame << " (t = " << t << "): "
           << bad_phase << " phase and " << bad_conc << " concentration values";
    record_failure(run_folder, reason.str());
    return false;
}

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

// The computation: the time stepping of u (the initial condition, set by start) with the model,
// and the frames written into run_folder. Returns the exit code of the program.
template< class ModelT >
int run(const Parameters& parameters, const ModelT& model, Vector& u, const RunFolder& run_folder)
{
    const Grid grid = parameters.grid();

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

    // Time of the previous attempted step (NaN: none yet) and the number of attempts in a row at it.
    Real last_attempt_time = std::numeric_limits<Real>::quiet_NaN();
    int attempts_at_same_time = 0;

    auto time_stepping = [ =, &rhs_evaluations, &last_attempt_time, &attempts_at_same_time ]
                         ( const Real& t, const Real& tau, const VectorView& u, VectorView& fu )
    {
        // The first stage of every attempt is evaluated at the solver time itself (time + 0 * tau),
        // so a repeated t there means a repeated attempt (rejected, or a step too small to move the time).
        if (rhs_evaluations % Method::getStages() == 0) {
            attempts_at_same_time = (t == last_attempt_time) ? attempts_at_same_time + 1 : 1;
            last_attempt_time = t;
            if (attempts_at_same_time > max_attempts_at_same_time) {
                std::ostringstream reason;
                reason << "the solver is stuck: " << max_attempts_at_same_time
                       << " attempts in a row at t = " << t << " (tau = " << tau << ")";
                throw SolverStuck(reason.str());
            }
        }
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
    solver.setAdaptivity( parameters.adaptivity );
    solver.setTime( parameters.initial_time );

    // Frame numbers are absolute (time / timeStep), so a run continued from a file with a later
    // initial_time carries on with the numbering of the previous run.
    const Index first_frame = static_cast<Index>(round(parameters.initial_time / parameters.timeStep));

    // Writes frame number `frame` (at time t). If that fails (full disk, quota, missing folder;
    // write_frame prints the details), records the failure and returns false: the run then stops
    // instead of going on without its frames.
    auto write_frame_or_fail = [ & ]( Index frame, Real t ) {
        if (write_frame(run_folder.calculations(), frame, u, grid))
            return true;
        std::ostringstream reason;
        reason << "unable to write frame " << frame << " (t = " << t << ")";
        record_failure(run_folder, reason.str());
        return false;
    };

    if (!write_frame_or_fail(first_frame, parameters.initial_time))
        return EXIT_FAILURE;
    // The initial condition too: a NaN in it (e.g. from a file) would otherwise only show up later.
    if (!check_frame_finite(u, grid, first_frame, parameters.initial_time, run_folder))
        return EXIT_RUN_FAILED;

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
        try {
            solver.solve( u, time_stepping );
        }
        catch (const SolverStuck& stuck) {
            // u is the last accepted state (TNL changes u only when it accepts a step, and the exception
            // comes before any stage is computed), somewhere between the previous frame and this one.
            const Index N = grid.num_points();
            std::ostringstream reason;
            reason << stuck.what() << "; NaN or inf values: " << count_non_finite(u, 0, N) << " phase, "
                   << count_non_finite(u, N, 2 * N) << " concentration";
            if (write_state(run_folder.stuck_state(), u, grid))
                reason << "; state saved in " << run_folder.stuck_state().string();
            record_failure(run_folder, reason.str());
            return EXIT_RUN_FAILED;  // the frames up to the previous one are plotted
        }
        if (!write_frame_or_fail(first_frame + i, solver.getTime()))
            return EXIT_FAILURE;
        auto stepEnd = std::chrono::high_resolution_clock::now();

        // solve() resets the iteration counter, so getIterations() is the steps accepted in this frame.
        // Every attempt (accepted or rejected) evaluates the RHS once per stage of the method.
        const long long frame_accepted = solver.getIterations();
        const long long frame_attempts = (rhs_evaluations - rhs_evaluations_before) / Method::getStages();
        accepted_steps += frame_accepted;

        int elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(stepEnd - stepStart).count();
        print_progress(i, parameters, elapsed, frame_accepted, frame_attempts - frame_accepted);

        // The broken frame stays on disk, so it shows where the solution broke down. "Computation time"
        // is not written, so cluster.sh does not count the run as finished.
        if (!check_frame_finite(u, grid, first_frame + i, solver.getTime(), run_folder))
            return EXIT_RUN_FAILED;
    }

    const double computation_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - computation_start).count();

    if (rhs_evaluations % Method::getStages() != 0) {
        const std::string message = std::to_string(rhs_evaluations) + " RHS evaluations is not a multiple of the "
            + std::to_string(Method::getStages()) + " stages of the method, so the number of rejected steps is wrong.";
        warn(message);
        run_folder.append_warnings({ message });
    }
    const long long rejected_steps = rhs_evaluations / Method::getStages() - accepted_steps;
    std::cout << "Integration steps: " << accepted_steps << " (rejected: " << rejected_steps << ")" << std::endl;
    // Before the computation time: that line marks a finished run (cluster.sh), so everything else is written by then.
    run_folder.append_integration_steps(accepted_steps, rejected_steps);

    std::cout << "Computation time: " << format_computation_time(computation_seconds) << std::endl;
    run_folder.append_computation_time(computation_seconds);
    return EXIT_SUCCESS;
}

// Whether the parameters of a model contain ksi and the anisotropy of an anisotropic Allen-Cahn phase
// equation, which the initial conditions HyperbolicTangent and WulffShape use.
template< class Params, class = void >
struct has_phase_interface : std::false_type {};
template< class Params >
struct has_phase_interface< Params, std::void_t< decltype(&Params::ksi), decltype(&Params::anisotropy) > >
    : std::true_type {};

// Reads the parameters of the model ModelT from its config section (strictly, so that a wrong config
// stops before anything is written), creates run_folder and records the setup there, builds the model
// and the initial condition, and runs the computation.
template< class ModelT >
int start(const Parameters& parameters, const RunFolder& run_folder,
          const std::filesystem::path& input_config_path, std::time_t run_start)
{
    const typename ModelT::Params model_params = ModelT::Params::load(parameters.model_config);
    // The run settings the model cannot work with (Model::check_parameters, or the model's own).
    ModelT::check_parameters(parameters);

    // ksi and the anisotropy of the model, for the initial conditions HyperbolicTangent and WulffShape.
    // A model without them cannot use those two (checked here, before anything is written).
    Real ic_ksi = 0;
    Anisotropy ic_anisotropy{0, 0, 0};
    if constexpr (has_phase_interface< typename ModelT::Params >::value) {
        ic_ksi = model_params.ksi;
        ic_anisotropy = model_params.anisotropy;
    }
    else if (!parameters.init_cond_from_file
             && (parameters.init_condition == ICType::HyperbolicTangent
                 || parameters.init_condition == ICType::WulffShape))
        throw std::runtime_error("The initial condition " + ICType_to_string(parameters.init_condition)
                                 + " needs ksi and the anisotropy of the model, which the model "
                                 + ModelT::name + " does not have.");

    if (!run_folder.create())
        return EXIT_FAILURE;

    // From here on the run folder belongs to this run, so an unexpected exception (a file that cannot be
    // written, a CUDA error, TNL's "time step must be positive" after an overflow, ...) is recorded in it
    // like the failures run detects itself. Before, an exception only means a wrong config: main reports it.
    bool computing = false;  // the time loop has started, so frames are written
    try {
        run_folder.save_setup(parameters, model_params.entries(), input_config_path);
        run_folder.write_run_start_info(run_start);
        // The warnings so far, mostly about the config (read before the run folder existed).
        run_folder.append_warnings(warnings());

        const ModelT model(parameters, model_params);
        Vector u(model.getDegreesOfFreedom());
        if (!set_initial_condition(u, parameters, ic_ksi, ic_anisotropy,
                                   model_params.initial_concentrations(), run_folder.calculations()))
            return EXIT_FAILURE;

        computing = true;
        return run(parameters, model, u, run_folder);
    }
    catch (const std::exception& e) {
        record_failure(run_folder, std::string("unexpected error: ") + e.what());
        // Once the time loop runs, the frames written so far are plotted, as after NaN or a stuck solver.
        return computing ? EXIT_RUN_FAILED : EXIT_FAILURE;
    }
}

int main(int argc, char** argv) {
    const std::time_t run_start = std::time(nullptr);

    // An uncaught exception would end the program with std::terminate (abort, exit code 134); instead its
    // message is printed and the program exits with EXIT_FAILURE. These are mostly config errors, found
    // before the run folder is touched; later exceptions are also recorded in the run folder by start.
    try {
        // usage: Phase-field [result_folder [config_file]]
        std::filesystem::path result_path = "Results";
        if (argc >= 2)
            result_path = argv[1];
        std::filesystem::path input_config_path = "config/config.json";
        if (argc >= 3)
            input_config_path = argv[2];

        const Parameters parameters = Parameters::load(input_config_path);
        const RunFolder run_folder(result_path);

        // The model is chosen by its name in the config (model.type). Each model is a separate
        // instantiation of start and run, compiled with its own GPU kernels; a new model is added here.
        if (parameters.model_type == AllenCahn::name)
            return start<AllenCahn>(parameters, run_folder, input_config_path, run_start);
        if (parameters.model_type == ZrNbChoud2013::name)
            return start<ZrNbChoud2013>(parameters, run_folder, input_config_path, run_start);
        if (parameters.model_type == ZrHMa2002::name)
            return start<ZrHMa2002>(parameters, run_folder, input_config_path, run_start);
        throw std::runtime_error("Unknown model in the config: model.type = \"" + parameters.model_type + "\" (known: "
                                 + AllenCahn::name + ", " + ZrNbChoud2013::name + ", " + ZrHMa2002::name + ").");
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }
    catch (...) {
        std::cerr << "Error: unknown exception." << std::endl;
        return EXIT_FAILURE;
    }
}