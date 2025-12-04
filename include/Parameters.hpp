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
    // initial condition parameters
    ICType init_condition;
    bool init_cond_from_file;
    std::string init_cond_file_path;
    double init_cond_radius;
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
    double m;
    double theta_0;

    // --- metody ---
    Parameters() = default;

    // načtení z JSON souboru
    static Parameters load(const std::filesystem::path& filename);

    // uložení do textového souboru (pro info/parameters.txt)
    void save_human_readable(const std::filesystem::path& filename) const;

    // uložení do textového souboru ve formátu vhodném pro zkopírování do LateXu
    void save_for_latex(const std::filesystem::path& filename) const;

    // uloží přesnou kopii původního JSON
    void save_copy_of_config(const std::filesystem::path& original_path,
                             const std::filesystem::path& copy_path) const;

};