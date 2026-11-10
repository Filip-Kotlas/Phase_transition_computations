#pragma once

#include <string>
#include <json.hpp>
#include <filesystem>
#include "constants.hpp"
#include "enums.hpp"

class Parameters {
public:
    // solver parameters
    double initial_time;
    double final_time;
    Domain domain;
    int sizeX;
    int sizeY;
    int frame_num;
    double timeStep;
    double integrationTimeStep;
    OutputFormat output_format;
    // initial condition parameters
    ICType init_condition;
    bool init_cond_from_file;
    std::string init_cond_file_path;
    double init_cond_radius;
    double init_cond_slope_width;
    bool init_cond_scaling;
    // boundary conditions
    BCType bc_phase_x;
    BCType bc_phase_y;
    BCType bc_conc_x;
    BCType bc_conc_y;
    // model parameters
    double alpha;
    double beta;
    double par_a;
    double par_b;
    double par_d;
    double T;
    double ksi;
    // force term parameters
    FTType force_term_type;
    double force_term_size;
    // anisotropy parameters
    double A;
    int m;
    double theta_0;

    // --- metody ---
    Parameters() = default;

    // load from a JSON file
    static Parameters load(const std::filesystem::path& filename);

    // save to a text file (for info/setup/parameters.txt)
    void save_human_readable(const std::filesystem::path& filename) const;

    // save to a text file in a format suitable for copying into LaTeX
    void save_for_latex(const std::filesystem::path& filename) const;

    // save an exact copy of the original JSON
    void save_copy_of_config(const std::filesystem::path& original_path,
                             const std::filesystem::path& copy_path) const;

};