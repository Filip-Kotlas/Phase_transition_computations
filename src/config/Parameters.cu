#include "config/Parameters.hpp"
#include "config/ConfigSection.hpp"
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

    // The sections solver, init_cond and bound_cond are read strictly, like the model's section
    // (ConfigSection): every key is required, so a typo is an error instead of a silent default.
    for (const char* section : { "solver", "init_cond", "bound_cond" }) {
        if (!config.contains(section))
            throw std::runtime_error(std::string("The config is missing the section ") + section + ".");
    }

    // solver parameters
    ConfigSection solver(config.at("solver"), "solver");
    p.initial_time  = solver.get<Real>("initial_time");
    p.final_time    = solver.get<Real>("final_time");
    p.sizeX         = solver.get<int>("sizeX");
    p.sizeY         = solver.get<int>("sizeY");
    p.frame_num     = solver.get<int>("frame_num");
    if (p.sizeX < 3 || p.sizeY < 3)
        throw std::runtime_error("The grid needs at least 3 points in each direction (solver.sizeX, solver.sizeY).");
    if (p.frame_num < 1)
        throw std::runtime_error("solver.frame_num has to be at least 1.");
    if (!(p.final_time > p.initial_time))
        throw std::runtime_error("solver.final_time has to be larger than solver.initial_time.");
    p.timeStep      = (p.final_time - p.initial_time) / p.frame_num;

    // domain
    ConfigSection domain(solver.get<json>("domain"), "solver.domain");
    p.domain.x_left  = domain.get<Real>("x_left");
    p.domain.x_right = domain.get<Real>("x_right");
    p.domain.y_left  = domain.get<Real>("y_left");
    p.domain.y_right = domain.get<Real>("y_right");
    domain.check_all_read();
    if (!(p.domain.x_right > p.domain.x_left) || !(p.domain.y_right > p.domain.y_left))
        throw std::runtime_error("The domain has to have x_right > x_left and y_right > y_left (solver.domain).");

    // integration time step (integration_time_step is required even when it is not used, like every key)
    const Grid grid = p.grid();
    Real computed_integration_step = pow(std::min(grid.hx, grid.hy), 2) / 5.0;
    const bool custom_integration_step = solver.get<bool>("custom_integration_time_step");
    const Real integration_time_step = solver.get<Real>("integration_time_step");
    p.adaptivity = solver.get<Real>("adaptivity");
    solver.check_all_read();
    if (!(p.adaptivity >= 0))
        throw std::runtime_error("solver.adaptivity has to be 0 (no adaptivity) or positive.");
    if (custom_integration_step) {
        if (!(integration_time_step > 0))
            throw std::runtime_error("solver.integration_time_step has to be positive.");
        p.integrationTimeStep = integration_time_step;
    } else {
        p.integrationTimeStep = computed_integration_step;
    }

    // initial condition parameters
    ConfigSection initial_condition(config.at("init_cond"), "init_cond");
    p.init_cond_from_file = initial_condition.get<bool>("init_cond_from_file");
    p.init_cond_file_path = initial_condition.get<std::string>("init_cond_file_path");
    p.init_condition = string_to_ICType(initial_condition.get<std::string>("initial_condition"));
    Real smaller_side = std::min(p.domain.x_right - p.domain.x_left,
                                 p.domain.y_right - p.domain.y_left);
    // Read as double and multiplied in double, then rounded to Real once (as before), so the radius
    // does not change in the last bit.
    p.init_cond_radius = smaller_side * initial_condition.get<double>("radius_proportion");
    p.init_cond_slope_width = initial_condition.get<Real>("slope_width");
    p.init_cond_scaling = initial_condition.get<bool>("scaling");
    initial_condition.check_all_read();

    // boundary conditions
    ConfigSection boundary_conditions(config.at("bound_cond"), "bound_cond");
    p.bound_cond.phase_x = string_to_BCType(boundary_conditions.get<std::string>("x_phase"));
    p.bound_cond.phase_y = string_to_BCType(boundary_conditions.get<std::string>("y_phase"));
    p.bound_cond.conc_x  = string_to_BCType(boundary_conditions.get<std::string>("x_conc"));
    p.bound_cond.conc_y  = string_to_BCType(boundary_conditions.get<std::string>("y_conc"));
    boundary_conditions.check_all_read();
    if ( (p.bound_cond.phase_x != BCType::Neumann) && (p.bound_cond.phase_x != BCType::Dirichlet) )
        throw std::runtime_error("Unknown boundary condition for phase field in x direction: " + BCType_to_string(p.bound_cond.phase_x));
    if ( (p.bound_cond.phase_y != BCType::Neumann) && (p.bound_cond.phase_y != BCType::Dirichlet) )
        throw std::runtime_error("Unknown boundary condition for phase field in y direction: " + BCType_to_string(p.bound_cond.phase_y));
    if ( (p.bound_cond.conc_x != BCType::Neumann) && (p.bound_cond.conc_x != BCType::Dirichlet) )
        throw std::runtime_error("Unknown boundary condition for concentration in x direction: " + BCType_to_string(p.bound_cond.conc_x));
    if ( (p.bound_cond.conc_y != BCType::Neumann) && (p.bound_cond.conc_y != BCType::Dirichlet) )
        throw std::runtime_error("Unknown boundary condition for concentration in y direction: " + BCType_to_string(p.bound_cond.conc_y));

    // the model: only its name (model) and its section (parameters) are taken here, the model reads
    // the section itself (ModelT::Params::load), strictly like the sections above.
    if (!config.contains("model"))
        throw std::runtime_error("The config is missing the entry model (the name of the model).");
    if (!config.at("model").is_string())
        throw std::runtime_error("The config entry model has to be the name of the model, it is "
                                 + config.at("model").dump() + ".");
    p.model_type = config.at("model").get<std::string>();
    if (!config.contains("parameters"))
        throw std::runtime_error("The config is missing the section parameters (the physical parameters of the model "
                                 + p.model_type + ").");
    p.model_config = config.at("parameters");

    // Nothing else is allowed at the top level, e.g. a parameter that belongs into parameters.
    for (const auto& item : config.items()) {
        const std::string& key = item.key();
        if (key != "solver" && key != "init_cond" && key != "bound_cond" && key != "model" && key != "parameters")
            throw std::runtime_error("The config section " + key + " is unknown (the physical parameters"
                                     " of the model belong into parameters).");
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
    file << std::left << std::setw(max_label_length) << "Adaptivity:" << std::right << std::setw(max_value_length) << adaptivity << std::endl;
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
    for (const ParameterEntry& entry : model_text_entries())
        file << std::left << std::setw(max_label_length) << entry.label + ":" << std::right << std::setw(max_value_length) << entry.value << std::endl;
}

namespace {
// Text for LaTeX outside math mode: a bare _ is an error there.
std::string latex_text(const std::string& text)
{
    std::string result;
    for (const char ch : text) {
        if (ch == '_')
            result += '\\';
        result += ch;
    }
    return result;
}
}

std::vector<ParameterEntry> Parameters::model_text_entries() const
{
    std::vector<ParameterEntry> result;
    for (const auto& item : model_config.items()) {
        if (item.value().is_string())
            result.push_back({item.key(), latex_text(item.key()), latex_text(item.value().get<std::string>())});
    }
    return result;
}

std::filesystem::path Parameters::model_database() const
{
    if (model_config.contains("database") && model_config.at("database").is_string())
        return model_config.at("database").get<std::string>();
    return {};
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
    file << "adaptivity & " << adaptivity << " \\\\" << std::endl;
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
    for (const ParameterEntry& entry : model_text_entries())
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
        throw std::runtime_error("Unable to write the copy of the config: " + copy_path.string());
    }
    dst << src.rdbuf();
}