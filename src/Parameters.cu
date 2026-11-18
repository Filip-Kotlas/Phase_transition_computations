#include "Parameters.hpp"
#include "ConfigSection.hpp"
#include <fstream>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <cmath>
#include <algorithm>

using nlohmann::json;

Parameters Parameters::load(const std::filesystem::path& filename) {
    std::ifstream file(filename);
    if (!file) {
        throw std::runtime_error("Unable to open the file " + filename.string());
    }

    json config;
    file >> config;
    Parameters p;

    // solver parameters
    auto solver = config.at("solver");
    p.initial_time  = solver.value("initial_time", 0.0);
    p.final_time    = solver.value("final_time", 0.3);
    p.sizeX         = solver.value("sizeX", 200);
    p.sizeY         = solver.value("sizeY", 200);
    p.frame_num     = solver.value("frame_num", 100);
    p.timeStep      = (p.final_time - p.initial_time) / p.frame_num;

    // domain
    auto domain = solver.at("domain");
    p.domain.x_left  = domain.value("x_left", 0.0);
    p.domain.x_right = domain.value("x_right", 1.0);
    p.domain.y_left  = domain.value("y_left", 0.0);
    p.domain.y_right = domain.value("y_right", 1.0);
   
    // integration time step
    const Grid grid = p.grid();
    Real computed_integration_step = pow(std::min(grid.hx, grid.hy), 2) / 5.0;
    if (solver.value("custom_integration_time_step", false)) {
        p.integrationTimeStep = solver.value("integration_time_step", computed_integration_step);
    } else {
        p.integrationTimeStep = computed_integration_step;
    }

    // initial condition parameters
    auto initial_condition = config.at("init_cond");
    p.init_cond_from_file = initial_condition.value("init_cond_from_file", false);
    p.init_cond_file_path = initial_condition.value("init_cond_file_path", "");
    p.init_condition = string_to_ICType(initial_condition.value("initial_condition", "constant_circle"));
    Real smaller_side = std::min(p.domain.x_right - p.domain.x_left,
                                 p.domain.y_right - p.domain.y_left);
    p.init_cond_radius = smaller_side * initial_condition.value("radius_proportion", 0.5);
    p.init_cond_slope_width = initial_condition.value("slope_width", 0.05);
    p.init_cond_scaling = initial_condition.value("scaling", false);

    // boundary conditions
    auto boundary_conditions = config.at("bound_cond");
    p.bound_cond.phase_x = string_to_BCType(boundary_conditions.value("x_phase", "neumann"));
    p.bound_cond.phase_y = string_to_BCType(boundary_conditions.value("y_phase", "neumann"));
    p.bound_cond.conc_x  = string_to_BCType(boundary_conditions.value("x_conc", "neumann"));
    p.bound_cond.conc_y  = string_to_BCType(boundary_conditions.value("y_conc", "neumann"));
    if ( (p.bound_cond.phase_x != BCType::Neumann) && (p.bound_cond.phase_x != BCType::Dirichlet) )
        throw std::runtime_error("Unknown boundary condition for phase field in x direction: " + BCType_to_string(p.bound_cond.phase_x));
    if ( (p.bound_cond.phase_y != BCType::Neumann) && (p.bound_cond.phase_y != BCType::Dirichlet) )
        throw std::runtime_error("Unknown boundary condition for phase field in y direction: " + BCType_to_string(p.bound_cond.phase_y));
    if ( (p.bound_cond.conc_x != BCType::Neumann) && (p.bound_cond.conc_x != BCType::Dirichlet) )
        throw std::runtime_error("Unknown boundary condition for concentration in x direction: " + BCType_to_string(p.bound_cond.conc_x));
    if ( (p.bound_cond.conc_y != BCType::Neumann) && (p.bound_cond.conc_y != BCType::Dirichlet) )
        throw std::runtime_error("Unknown boundary condition for concentration in y direction: " + BCType_to_string(p.bound_cond.conc_y));

    // The old layout had the model's parameters in the top-level sections problem, force_term and anisotropy.
    if (!config.contains("model") && (config.contains("problem") || config.contains("force_term")))
        throw std::runtime_error("The config " + filename.string() + " has the old layout (the sections problem,"
                                 " force_term, anisotropy); convert it with python3 python/migrate_config.py "
                                 + filename.string());

    // the model: only its type and its section are taken here, the model reads the section itself
    // (ModelT::Params::load). The section "model" may also hold the sections of other models.
    if (!config.contains("model"))
        throw std::runtime_error("The config is missing the section model.");
    ConfigSection model(config.at("model"), "model");
    p.model_type = model.get<std::string>("type");
    if (!model.contains(p.model_type))
        throw std::runtime_error("The config is missing model." + p.model_type
                                 + ", the section with the parameters of the model " + p.model_type + ".");
    p.model_config = config.at("model").at(p.model_type);

    // Anything else at the top level is not read, e.g. a parameter that belongs into model.<type>.
    for (const auto& item : config.items()) {
        if (item.key() != "solver" && item.key() != "init_cond" && item.key() != "bound_cond" && item.key() != "model")
            std::cerr << "Warning: the config section " << item.key() << " is not used (the physical parameters"
                         " of the model belong into model." << p.model_type << ")." << std::endl;
    }

    return p;
}

void Parameters::save_human_readable(const std::filesystem::path& filename,
                                     const std::vector<ParameterEntry>& model_entries) const {
    std::ofstream file(filename);
    if (!file) {
        throw std::runtime_error("Unable to open the file " + filename.string());
    }

    int max_label_length = 30;
    int max_value_length = 28;

    file << std::left << std::setw(max_label_length) << "Initial time:" << std::right << std::setw(max_value_length) << initial_time << std::endl;
    file << std::left << std::setw(max_label_length) << "Final time:"   << std::right << std::setw(max_value_length) << final_time << std::endl;
    
    std::ostringstream oss;
    oss << "[(" << std::fixed << std::setprecision(2) 
        << domain.x_left << ", " << domain.x_right << ")"
        << "(" << domain.y_left << ", " << domain.y_right << ")]";
    file << std::left << std::setw(max_label_length) << "Domain:" << std::right << std::setw(max_value_length) << oss.str() << std::endl;
    file << std::left << std::setw(max_label_length) << "SizeX:"  << std::right << std::setw(max_value_length) << sizeX << std::endl;
    file << std::left << std::setw(max_label_length) << "SizeY:"  << std::right << std::setw(max_value_length) << sizeY << std::endl;
    file << std::left << std::setw(max_label_length) << "Time step:" << std::right << std::setw(max_value_length) << timeStep << std::endl;
    file << std::left << std::setw(max_label_length) << "Integration time step:" << std::right << std::setw(max_value_length) << integrationTimeStep << std::endl;
    file << std::endl;
    file << std::left << std::setw(max_label_length) << "Initial condition type:"   << std::right << std::setw(max_value_length) << ICType_to_string(init_condition) << std::endl;
    file << std::left << std::setw(max_label_length) << "Initial condition radius:"   << std::right << std::setw(max_value_length) << init_cond_radius << std::endl;
    file << std::left << std::setw(max_label_length) << "Initial condition slope width:"   << std::right << std::setw(max_value_length) << init_cond_slope_width << std::endl;
    file << std::endl;
    file << std::left << std::setw(max_label_length) << "Boundary condition phase x:"   << std::right << std::setw(max_value_length) << BCType_to_string(bound_cond.phase_x) << std::endl;
    file << std::left << std::setw(max_label_length) << "Boundary condition phase y:"   << std::right << std::setw(max_value_length) << BCType_to_string(bound_cond.phase_y) << std::endl;
    file << std::left << std::setw(max_label_length) << "Boundary condition conc x:"    << std::right << std::setw(max_value_length) << BCType_to_string(bound_cond.conc_x) << std::endl;
    file << std::left << std::setw(max_label_length) << "Boundary condition conc y:"    << std::right << std::setw(max_value_length) << BCType_to_string(bound_cond.conc_y) << std::endl;
    file << std::endl;
    file << std::left << std::setw(max_label_length) << "Model:"   << std::right << std::setw(max_value_length) << model_type << std::endl;
    for (const ParameterEntry& entry : model_entries)
        file << std::left << std::setw(max_label_length) << entry.label + ":" << std::right << std::setw(max_value_length) << entry.value << std::endl;
}

void Parameters::save_for_latex(const std::filesystem::path& filename,
                                const std::vector<ParameterEntry>& model_entries) const {
    std::ofstream file(filename);
    if (!file) {
        throw std::runtime_error("Unable to open the file " + filename.string());
    }

    file << "\\textbf{Domain:}" << std::endl << std::endl;
    file << "\\begin{tabular}{ll}" << std::endl;
    file << "\\(\\Omega\\) & \\((" << std::defaultfloat 
    << domain.x_left << ", " << domain.x_right << ") \\times "
    << "(" << domain.y_left << ", " << domain.y_right << ")\\) \\\\" << std::endl;
    file << "\\(N_x\\) & " << sizeX << " \\\\" << std::endl;
    file << "\\(N_y\\) & " << sizeY << " \\\\" << std::endl;
    file << "\\end{tabular}" << std::endl << std::endl;

    file << "\\textbf{Time parameters:}" << std::endl << std::endl;
    file << "\\begin{tabular}{ll}" << std::endl;
    file << "\\(t_{0}\\) & " << initial_time << " \\\\" << std::endl;
    file << "\\(t_{max}\\) & " << final_time << " \\\\" << std::endl;
    file << "\\(\\tau\\) & " << integrationTimeStep << " \\\\" << std::endl;
    file << "\\end{tabular}" << std::endl << std::endl;

    file << "\\textbf{Initial condition:}" << std::endl << std::endl;
    file << "\\begin{tabular}{ll}" << std::endl;
    file << "type & " << ICType_to_string(init_condition) << " \\\\" << std::endl;
    file << "radius & " << init_cond_radius << " \\\\" << std::endl;
    file << "transition width & " << init_cond_slope_width << " \\\\" << std::endl;
    file << "\\end{tabular}" << std::endl << std::endl;

    file << "\\textbf{Boundary conditions:}" << std::endl << std::endl;
    file << "\\begin{tabular}{ll}" << std::endl;
    if (bound_cond.phase_x == BCType::Dirichlet && bound_cond.phase_y == BCType::Dirichlet &&
        bound_cond.conc_x == BCType::Dirichlet && bound_cond.conc_y == BCType::Dirichlet) {
        file << "everywhere & Dirichlet \\\\" << std::endl << std::endl;
    }
    else if (bound_cond.phase_x == BCType::Neumann && bound_cond.phase_y == BCType::Neumann &&
             bound_cond.conc_x == BCType::Neumann && bound_cond.conc_y == BCType::Neumann) {
        file << "everywhere & Neumann \\\\" << std::endl << std::endl;
    }
    else {            
        file << "phase field in the \\(x\\) direction & " << BCType_to_string(bound_cond.phase_x) << " \\\\" << std::endl;
        file << "phase field in the \\(y\\) direction & " << BCType_to_string(bound_cond.phase_y) << " \\\\" << std::endl;
        file << "concentration in the \\(x\\) direction & " << BCType_to_string(bound_cond.conc_x) << " \\\\" << std::endl;
        file << "concentration in the \\(y\\) direction & " << BCType_to_string(bound_cond.conc_y) << " \\\\" << std::endl;
    }
    file << "\\end{tabular}" << std::endl << std::endl;

    std::string latex_model_type;  // "allen_cahn" -> "allen\_cahn", a bare _ is an error outside math in LaTeX
    for (const char ch : model_type) {
        if (ch == '_')
            latex_model_type += '\\';
        latex_model_type += ch;
    }
    file << "\\textbf{Model parameters (" << latex_model_type << "):}" << std::endl << std::endl;
    file << "\\begin{tabular}{ll}" << std::endl;
    for (const ParameterEntry& entry : model_entries)
        file << entry.latex << " & " << entry.value << " \\\\" << std::endl;
    file << "\\end{tabular}" << std::endl;
}

void Parameters::save_copy_of_config(const std::filesystem::path& original_path,
                                     const std::filesystem::path& copy_path) const {
    std::ifstream src(original_path, std::ios::binary);
    if (!src) {
        throw std::runtime_error("Unable to open original config: " + original_path.string());
    }
    std::ofstream dst(copy_path, std::ios::binary);
    if (!dst) {
        throw std::runtime_error("Unable to open original config: " + copy_path.string());
    }
    dst << src.rdbuf();
}