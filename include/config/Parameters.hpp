#pragma once

#include <cmath>
#include <string>
#include <vector>
#include <sstream>
#include "external/json.hpp"
#include <filesystem>
#include "core/types.hpp"
#include "config/enums.hpp"
#include "core/Grid.hpp"

// Plain data that the GPU code needs (the models, InitialCondition): numbers and enums only, so that
// copying them to the GPU by value is well defined. Parameters itself, which also holds strings and
// the settings of the run, stays on the host.

// The anisotropy psi(theta) = 1 + A cos(m (theta - theta_0)) of the interface energy (a member of the
// Params of AllenCahnModel and ZirconiumModel, also used by the initial condition WulffShape).
struct Anisotropy {
    Real A;
    int m;
    Real theta_0;

    __cuda_callable__
    Real psi(const Real theta) const { return 1 + A*cos(m*(theta - theta_0)); }

    // d psi / d theta
    __cuda_callable__
    Real der_psi(const Real theta) const { return -A*m*sin(m*(theta - theta_0)); }
};

// The boundary conditions, common to all models (the config section bound_cond, read by
// Parameters::load into Parameters::bound_cond; the model and the initial condition copy them).
struct BoundaryConditions {
    BCType phase_x;
    BCType phase_y;
    BCType conc_x;
    BCType conc_y;
};

// One parameter of a model, as written to parameters.txt (label: value) and parameters.tex
// (latex & value). Each model lists its own parameters (ModelT::Params::entries()).
struct ParameterEntry {
    std::string label;   // e.g. "T"
    std::string latex;   // e.g. \(T\)
    std::string value;   // already formatted, see to_text

    template< class T >
    static std::string to_text(const T& value)
    {
        std::ostringstream text;
        text << value;
        return text.str();
    }
};

// The configuration of a run (the config file, see Parameters::load): the solver, the initial
// condition, the boundary conditions and which model to solve. The physical parameters of the model
// are read by the model itself (ModelT::Params::load) from model_config. Host only.
struct Parameters {
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
    // boundary conditions (plain data, copied to the GPU)
    BoundaryConditions bound_cond;
    // the model (the config section "model"): its name (model.type, e.g. "zirconium") chooses the
    // model class in main.cu, its section model.<type> holds the physical parameters of that model
    std::string model_type;
    nlohmann::json model_config;

    Parameters() = default;

    // The grid given by sizeX, sizeY and domain.
    Grid grid() const { return Grid(sizeX, sizeY, domain); }

    // load from a JSON file
    static Parameters load(const std::filesystem::path& filename);

    // save to a text file (for info/setup/parameters.txt), with the parameters of the model
    void save_human_readable(const std::filesystem::path& filename,
                             const std::vector<ParameterEntry>& model_entries) const;

    // save to a text file in a format suitable for copying into LaTeX, with the parameters of the model
    void save_for_latex(const std::filesystem::path& filename,
                        const std::vector<ParameterEntry>& model_entries) const;

    // save an exact copy of the original JSON
    void save_copy_of_config(const std::filesystem::path& original_path,
                             const std::filesystem::path& copy_path) const;

};