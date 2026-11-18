#include "RunFolder.hpp"

#include <cmath>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <unistd.h>
#include <TNL/Backend/Functions.h>
#include <TNL/Backend/DeviceInfo.h>

namespace {

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

}  // namespace

RunFolder::RunFolder(const std::filesystem::path& root)
:  calculations_dir(root / "calculations"),
   setup_dir(root / "info" / "setup"),
   parameters_file(setup_dir / "parameters.txt")
{}

bool RunFolder::create() const
{
    try {
        std::filesystem::create_directories(setup_dir);
        std::filesystem::create_directories(calculations_dir);
    }
    catch (const std::filesystem::filesystem_error& e) {
        std::cerr << "Error creating folders: " << e.what() << std::endl;
        return false;
    }
    return true;
}

void RunFolder::save_setup(const Parameters& parameters, const std::vector<ParameterEntry>& model_entries,
                           const std::filesystem::path& input_config) const
{
    parameters.save_human_readable(parameters_file, model_entries);
    parameters.save_copy_of_config(input_config, setup_dir / "config.json");
    parameters.save_for_latex(setup_dir / "parameters.tex", model_entries);
}

void RunFolder::write_run_start_info(std::time_t run_start) const
{
    std::ofstream file(parameters_file, std::ios::app);
    file << "\nRun date: " << std::put_time(std::localtime(&run_start), "%d.%m.%Y %H:%M:%S") << '\n'
         << "Computer: " << host_name() << '\n'
         << "GPU: " << gpu_description() << '\n';
    // Checked after close(), which writes the buffer to the disk, so that an error there is reported too.
    file.close();
    if (!file)
        std::cerr << "Warning: unable to write the run date, computer and GPU to " << parameters_file << std::endl;
}

void RunFolder::append_integration_steps(long long accepted, long long rejected) const
{
    std::ofstream file(parameters_file, std::ios::app);
    file << "Integration steps: " << accepted << " (rejected: " << rejected << ")\n";
    file.close();
    if (!file)
        std::cerr << "Warning: unable to write the number of integration steps to " << parameters_file << std::endl;
}

void RunFolder::append_computation_time(double seconds) const
{
    std::ofstream file(parameters_file, std::ios::app);
    file << "Computation time: " << format_computation_time(seconds) << '\n';
    file.close();
    if (!file)
        std::cerr << "Warning: unable to write the computation time to " << parameters_file << std::endl;
}

std::string format_computation_time(double seconds)
{
    const long total_seconds = std::lround(seconds);
    std::ostringstream text;
    text << total_seconds / 3600 << "h "
         << std::setfill('0') << std::setw(2) << (total_seconds % 3600) / 60 << "m "
         << std::setw(2) << total_seconds % 60 << "s"
         << " (" << std::fixed << std::setprecision(1) << seconds << " s)";
    return text.str();
}
