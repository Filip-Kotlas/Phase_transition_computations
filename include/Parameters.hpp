#pragma once

#include <string>
#include <json.hpp>
#include <filesystem>
#include "types.hpp"
#include "constants.hpp"
#include "enums.hpp"

class Parameters {
public:
    // solver parameters
    Real initial_time;
    Real final_time;
    Domain domain;
    int sizeX;
    int sizeY;
    int frame_num;
    Real timeStep;
    Real integrationTimeStep;
    // initial condition parameters
    ICType init_condition;
    bool init_cond_from_file;
    std::string init_cond_file_path;
    Real init_cond_radius;
    Real init_cond_slope_width;
    bool init_cond_scaling;
    // boundary conditions
    BCType bc_phase_x;
    BCType bc_phase_y;
    BCType bc_conc_x;
    BCType bc_conc_y;
    // model parameters
    Real alpha;
    Real beta;
    Real par_a;
    Real par_b;
    Real par_d;
    Real T;
    Real ksi;
    // force term parameters
    FTType force_term_type;
    Real force_term_size;
    // anisotropy parameters
    Real A;
    int m;
    Real theta_0;

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