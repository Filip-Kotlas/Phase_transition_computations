#include "Parameters.hpp"
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
    double computed_integration_step = pow(std::min((p.domain.x_right - p.domain.x_left) / (p.sizeX - 1),
                                                    (p.domain.y_right - p.domain.y_left) / (p.sizeY - 1)),
                                           2) / 5.0;
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
    double smaller_side = std::min(p.domain.x_right - p.domain.x_left,
                                   p.domain.y_right - p.domain.y_left);
    p.init_cond_radius = smaller_side * initial_condition.value("radius_proportion", 0.5);
    p.init_cond_slope_width = initial_condition.value("slope_width", 0.05);
    p.init_cond_scaling = initial_condition.value("scaling", false);

    // boundary conditions
    auto boundary_conditions = config.at("bound_cond");
    p.bc_phase_x = string_to_BCType(boundary_conditions.value("x_phase", "neumann"));
    p.bc_phase_y = string_to_BCType(boundary_conditions.value("y_phase", "neumann"));
    p.bc_conc_x  = string_to_BCType(boundary_conditions.value("x_conc", "neumann"));
    p.bc_conc_y  = string_to_BCType(boundary_conditions.value("y_conc", "neumann"));
    if ( (p.bc_phase_x != BCType::Neumann) && (p.bc_phase_x != BCType::Dirichlet) )
        throw std::runtime_error("Unknown boundary condition for phase field in x direction: " + BCType_to_string(p.bc_phase_x));
    if ( (p.bc_phase_y != BCType::Neumann) && (p.bc_phase_y != BCType::Dirichlet) )
        throw std::runtime_error("Unknown boundary condition for phase field in y direction: " + BCType_to_string(p.bc_phase_y));
    if ( (p.bc_conc_x != BCType::Neumann) && (p.bc_conc_x != BCType::Dirichlet) )
        throw std::runtime_error("Unknown boundary condition for concentration in x direction: " + BCType_to_string(p.bc_conc_x));
    if ( (p.bc_conc_y != BCType::Neumann) && (p.bc_conc_y != BCType::Dirichlet) )
        throw std::runtime_error("Unknown boundary condition for concentration in y direction: " + BCType_to_string(p.bc_conc_y));

    // problem parameters
    auto problem = config.at("problem");
    p.alpha = problem.value("alpha", 1.0);
    p.beta  = problem.value("beta", 1.0);
    p.par_a = problem.value("a", 1.0);
    p.par_b = problem.value("b", 0.1);
    p.par_d = problem.value("d", 5e15);
    p.T = problem.value("T", 1200);
    p.ksi   = problem.value("ksi", 0.01);

    // force term parameters
    auto force_term = config.at("force_term");
    p.force_term_type = string_to_FTType(force_term.value("type", "zirconium"));
    p.force_term_size = force_term.value("size", 0.0);
    if (p.force_term_type != FTType::Zirconium && p.force_term_type != FTType::Constant && p.force_term_type != FTType::Radial) {
        throw std::runtime_error("Unknown force term type in config: " + FTType_to_string(p.force_term_type));
    }

    // anisotropy parameters
    auto anisotropy = config.at("anisotropy");
    p.A = anisotropy.value("A", 0.3);
    p.m = anisotropy.value("m", 2);
    p.theta_0 = anisotropy.value("theta_0", 0.0);
    if (p.m <= 1)
        throw std::runtime_error("The number m is lower than 2.");
    if (p.A > 1.0/(p.m*p.m - 1))
        throw std::runtime_error("Parameters A and m do not fulfill the convexity condition.");
    
    return p;
}

void Parameters::save_human_readable(const std::filesystem::path& filename) const {
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
    file << std::left << std::setw(max_label_length) << "Boundary condition phase x:"   << std::right << std::setw(max_value_length) << BCType_to_string(bc_phase_x) << std::endl;
    file << std::left << std::setw(max_label_length) << "Boundary condition phase y:"   << std::right << std::setw(max_value_length) << BCType_to_string(bc_phase_y) << std::endl;
    file << std::left << std::setw(max_label_length) << "Boundary condition conc x:"    << std::right << std::setw(max_value_length) << BCType_to_string(bc_conc_x) << std::endl;
    file << std::left << std::setw(max_label_length) << "Boundary condition conc y:"    << std::right << std::setw(max_value_length) << BCType_to_string(bc_conc_y) << std::endl;
    file << std::endl;
    file << std::left << std::setw(max_label_length) << "Force term type:"   << std::right << std::setw(max_value_length) << FTType_to_string(force_term_type) << std::endl;
    if (force_term_type == FTType::Constant || force_term_type == FTType::Radial)
        file << std::left << std::setw(max_label_length) << "Force term size:"   << std::right << std::setw(max_value_length) << force_term_size << std::endl;
    file << std::endl;
    file << std::left << std::setw(max_label_length) << "Alpha:" << std::right << std::setw(max_value_length) << alpha << std::endl;
    file << std::left << std::setw(max_label_length) << "Beta:"  << std::right << std::setw(max_value_length) << beta << std::endl;
    file << std::left << std::setw(max_label_length) << "a:" << std::right << std::setw(max_value_length) << par_a << std::endl;
    file << std::left << std::setw(max_label_length) << "b:" << std::right << std::setw(max_value_length) << par_b << std::endl;
    file << std::left << std::setw(max_label_length) << "d:" << std::right << std::setw(max_value_length) << par_d << std::endl;
    file << std::left << std::setw(max_label_length) << "T:" << std::right << std::setw(max_value_length) << T << std::endl;
    file << std::left << std::setw(max_label_length) << "Ksi:"   << std::right << std::setw(max_value_length) << ksi << std::endl;
    file << std::left << std::setw(max_label_length) << "A:"   << std::right << std::setw(max_value_length) << A << std::endl;
    file << std::left << std::setw(max_label_length) << "m:"   << std::right << std::setw(max_value_length) << m << std::endl;
    file << std::left << std::setw(max_label_length) << "theta_0:"   << std::right << std::setw(max_value_length) << theta_0 << std::endl;
}

void Parameters::save_for_latex(const std::filesystem::path& filename) const {
    std::ofstream file(filename);
    if (!file) {
        throw std::runtime_error("Unable to open the file " + filename.string());
    }

    file << "\\textbf{Oblast:}" << std::endl << std::endl;
    file << "\\begin{tabular}{ll}" << std::endl;
    file << "\\(\\Omega\\) & \\((" << std::defaultfloat 
    << domain.x_left << ", " << domain.x_right << ") \\times "
    << "(" << domain.y_left << ", " << domain.y_right << ")\\) \\\\" << std::endl;
    file << "\\(N_x\\) & " << sizeX << " \\\\" << std::endl;
    file << "\\(N_y\\) & " << sizeY << " \\\\" << std::endl;
    file << "\\end{tabular}" << std::endl << std::endl;

    file << "\\textbf{Časové parametry:}" << std::endl << std::endl;
    file << "\\begin{tabular}{ll}" << std::endl;
    file << "\\(t_{0}\\) & " << initial_time << " \\\\" << std::endl;
    file << "\\(t_{max}\\) & " << final_time << " \\\\" << std::endl;
    file << "\\(\\tau\\) & " << integrationTimeStep << " \\\\" << std::endl;
    file << "\\end{tabular}" << std::endl << std::endl;

    file << "\\textbf{Počáteční podmínka:}" << std::endl << std::endl;
    file << "\\begin{tabular}{ll}" << std::endl;
    file << "typ & " << ICType_to_string(init_condition) << " \\\\" << std::endl;
    file << "poloměr & " << init_cond_radius << " \\\\" << std::endl;
    file << "šířka přechodu & " << init_cond_slope_width << " \\\\" << std::endl;
    file << "\\end{tabular}" << std::endl << std::endl;

    file << "\\textbf{Okrajové podmínky:}" << std::endl << std::endl;
    file << "\\begin{tabular}{ll}" << std::endl;
    if (bc_phase_x == BCType::Dirichlet && bc_phase_y == BCType::Dirichlet &&
        bc_conc_x == BCType::Dirichlet && bc_conc_y == BCType::Dirichlet) {
        file << "všude & Dirichlet \\\\" << std::endl << std::endl;
    }
    else if (bc_phase_x == BCType::Neumann && bc_phase_y == BCType::Neumann &&
             bc_conc_x == BCType::Neumann && bc_conc_y == BCType::Neumann) {
        file << "všude & Neumann \\\\" << std::endl << std::endl;
    }
    else {            
        file << "fázové pole v \\(x\\)-ovém směru & " << BCType_to_string(bc_phase_x) << " \\\\" << std::endl;
        file << "fázové pole v \\(y\\)-ovém směru & " << BCType_to_string(bc_phase_y) << " \\\\" << std::endl;
        file << "koncentrace v \\(x\\)-ovém směru & " << BCType_to_string(bc_conc_x) << " \\\\" << std::endl;
        file << "koncentrace v \\(y\\)-ovém směru & " << BCType_to_string(bc_conc_y) << " \\\\" << std::endl;
    }
    file << "\\end{tabular}" << std::endl << std::endl;

    file << "\\textbf{Parametry modelu:}" << std::endl << std::endl;
    file << "\\begin{tabular}{ll}" << std::endl;
    file << "\\(\\alpha\\) & " << alpha << " \\\\" << std::endl;
    file << "\\(\\beta\\) & " << beta << " \\\\" << std::endl;
    file << "\\(a\\) & " << par_a << " \\\\" << std::endl;
    file << "\\(b\\) & " << par_b << " \\\\" << std::endl;
    file << "\\(d\\) & " << par_d << " \\\\" << std::endl;
    file << "\\(T\\) & " << T << " \\\\" << std::endl;
    file << "\\(\\xi\\) & " << ksi << " \\\\" << std::endl;
    file << "\\end{tabular}" << std::endl << std::endl;

    file << "\\textbf{Parametry anizotropie:}" << std::endl << std::endl;
    file << "\\begin{tabular}{ll}" << std::endl;
    file << "\\(A\\) & " << A << " \\\\" << std::endl;
    file << "\\(m\\) & " << m << " \\\\" << std::endl;
    file << "\\(\\theta_0\\) & " << theta_0 << " \\\\" << std::endl;
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