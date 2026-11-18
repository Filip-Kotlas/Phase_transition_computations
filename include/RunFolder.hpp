#pragma once

// The output folder of one run, results/<FOLDER>/ (host code). It knows the layout of the folder
// and writes what the binary records about the run into info/setup/. The frames themselves are
// written by write_frame (ResultIO.hpp) into calculations(). The Python counterpart, which reads
// the folder, is Dataset in python/dataset.py; the Python scripts add info/plots and info/boundary.

#include <ctime>
#include <filesystem>
#include <string>
#include <vector>

#include "Parameters.hpp"

class RunFolder {
public:
    // root is results/<FOLDER> (the first argument of the binary). Only computes the paths,
    // touches nothing on disk.
    explicit RunFolder(const std::filesystem::path& root);

    // Creates info/setup/ and calculations/. Prints the error and returns false on failure.
    bool create() const;

    // The folder with the frames Result-NNNNN.npy; init_cond_file_path is relative to it too.
    const std::filesystem::path& calculations() const { return calculations_dir; }

    // Writes parameters.txt, a copy of the config (config.json) and parameters.tex into info/setup/;
    // model_entries are the parameters of the model (ModelT::Params::entries()).
    // Throws std::runtime_error when a file cannot be opened (from the Parameters::save_* methods).
    void save_setup(const Parameters& parameters, const std::vector<ParameterEntry>& model_entries,
                    const std::filesystem::path& input_config) const;

    // Appends the run date (run_start), the computer and the GPU to parameters.txt. Called at the
    // start, so that a run that crashes or is killed still records them. Only warns on failure.
    void write_run_start_info(std::time_t run_start) const;

    // Appends "Integration steps: N (rejected: M)" to parameters.txt: the steps of the adaptive
    // Kutta-Merson method accepted over the whole run, and the attempts it threw away because the error
    // estimate was too large. Called before append_computation_time. Only warns on failure.
    void append_integration_steps(long long accepted, long long rejected) const;

    // Appends "Computation time: ..." to parameters.txt ("make run" appends the git commit after it).
    // Only warns on failure.
    void append_computation_time(double seconds) const;

private:
    std::filesystem::path calculations_dir;  // results/<FOLDER>/calculations
    std::filesystem::path setup_dir;         // results/<FOLDER>/info/setup
    std::filesystem::path parameters_file;   // results/<FOLDER>/info/setup/parameters.txt
};

// Wall time in the form "0h 01m 23s (83.4 s)", as written to parameters.txt.
std::string format_computation_time(double seconds);
