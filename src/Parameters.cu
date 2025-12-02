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
    std::string ic_str = initial_condition.value("initial_condition", "hyperbolic_tangent");
    if (ic_str == "hyperbolic_tangent")         p.init_condition = ICType::HyperbolicTangent;
    else if (ic_str == "linear_by_parts")       p.init_condition = ICType::LinearByParts;
    else if (ic_str == "circle")                p.init_condition = ICType::ConstantCircle;
    else if (ic_str == "half")                  p.init_condition = ICType::ConstantHalves;
    else if (ic_str == "stripe")                p.init_condition = ICType::Stripe;
    else if (ic_str == "two_bumps")             p.init_condition = ICType::TwoBumps;
    else if (ic_str == "three_bumps")           p.init_condition = ICType::ThreeBumps;
    else if (ic_str == "star")                  p.init_condition = ICType::Star;
    else if (ic_str == "perpendicular_stripes") p.init_condition = ICType::PerpendicularStripes;
    else if (ic_str == "box")                   p.init_condition = ICType::Box;
    else if (ic_str == "random_bumps")          p.init_condition = ICType::RandomBumps;
    else if (ic_str == "wulff_shape")           p.init_condition = ICType::WulffShape;
    else throw std::runtime_error("Unknown initial_condition in config: " + ic_str);
    double smaller_side = std::min(p.domain.x_right - p.domain.x_left,
                                   p.domain.y_right - p.domain.y_left);
    p.init_cond_radius = smaller_side * initial_condition.value("radius_proportion", 0.5);

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

    file << std::left << std::setw(24) << "Initial time:" << std::right << std::setw(28) << initial_time << std::endl;
    file << std::left << std::setw(24) << "Final time:"   << std::right << std::setw(28) << final_time << std::endl;
    
    std::ostringstream oss;
    oss << "[(" << std::fixed << std::setprecision(2) 
        << domain.x_left << ", " << domain.x_right << ")"
        << "(" << domain.y_left << ", " << domain.y_right << ")]";
    file << std::left << std::setw(24) << "Domain:" << std::right << std::setw(28) << oss.str() << std::endl;
    file << std::left << std::setw(24) << "SizeX:"  << std::right << std::setw(28) << sizeX << std::endl;
    file << std::left << std::setw(24) << "SizeY:"  << std::right << std::setw(28) << sizeY << std::endl;
    file << std::left << std::setw(24) << "Time step:" << std::right << std::setw(28) << timeStep << std::endl;
    file << std::left << std::setw(24) << "Integration time step:" << std::right << std::setw(28) << integrationTimeStep << std::endl;
    file << std::endl;
    file << std::left << std::setw(24) << "Initial condition radius:"   << std::right << std::setw(28) << init_cond_radius << std::endl;
    file << std::endl;
    file << std::left << std::setw(24) << "Boundary condition phase x:"   << std::right << std::setw(28) << BCType_to_string(bc_phase_x) << std::endl;
    file << std::left << std::setw(24) << "Boundary condition phase y:"   << std::right << std::setw(28) << BCType_to_string(bc_phase_y) << std::endl;
    file << std::left << std::setw(24) << "Boundary condition conc x:"    << std::right << std::setw(28) << BCType_to_string(bc_conc_x) << std::endl;
    file << std::left << std::setw(24) << "Boundary condition conc y:"    << std::right << std::setw(28) << BCType_to_string(bc_conc_y) << std::endl;
    file << std::endl;
    file << std::left << std::setw(24) << "Force term type:"   << std::right << std::setw(28) << FTType_to_string(force_term_type) << std::endl;
    if (force_term_type == FTType::Constant || force_term_type == FTType::Radial)
        file << std::left << std::setw(24) << "Force term size:"   << std::right << std::setw(28) << force_term_size << std::endl;
    file << std::endl;
    file << std::left << std::setw(24) << "Alpha:" << std::right << std::setw(28) << alpha << std::endl;
    file << std::left << std::setw(24) << "Beta:"  << std::right << std::setw(28) << beta << std::endl;
    file << std::left << std::setw(24) << "a:" << std::right << std::setw(28) << par_a << std::endl;
    file << std::left << std::setw(24) << "b:" << std::right << std::setw(28) << par_b << std::endl;
    file << std::left << std::setw(24) << "d:" << std::right << std::setw(28) << par_d << std::endl;
    file << std::left << std::setw(24) << "T:" << std::right << std::setw(28) << T << std::endl;
    file << std::left << std::setw(24) << "Ksi:"   << std::right << std::setw(28) << ksi << std::endl;
    file << std::left << std::setw(24) << "A:"   << std::right << std::setw(28) << A << std::endl;
    file << std::left << std::setw(24) << "m:"   << std::right << std::setw(28) << m << std::endl;
    file << std::left << std::setw(24) << "theta_0:"   << std::right << std::setw(28) << theta_0 << std::endl;
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
    file << "poloměr počáteční podmínky: & " << init_cond_radius << " \\\\" << std::endl;
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