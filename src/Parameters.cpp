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

    json j;
    file >> j;

    Parameters p;

    auto solver = j.at("solver");
    p.type          = solver.value("type", "Runge-Kutta");
    p.initial_time  = solver.value("initial_time", 0.0);
    p.final_time    = solver.value("final_time", 0.3);
    p.sizeX         = solver.value("sizeX", 200);
    p.sizeY         = solver.value("sizeY", 200);
    p.frame_num     = solver.value("frame_num", 100);
    p.timeStep      = (p.final_time - p.initial_time) / p.frame_num;

    auto domain = solver.at("domain");
    p.domain.x_left  = domain.value("x_left", -1.0);
    p.domain.x_right = domain.value("x_right", 1.0);
    p.domain.y_left  = domain.value("y_left", -1.0);
    p.domain.y_right = domain.value("y_right", 1.0);

    int model_value = solver.value("model", 1);
    if (!(model_value >= 1 && model_value <= 4)) {
        throw std::runtime_error("Neplatná hodnota 'model' v JSON.");
    }
    p.model = static_cast<MODEL>(model_value);

    double computed_dt = pow(std::min((p.domain.x_right - p.domain.x_left) / (p.sizeX - 1),
                                      (p.domain.y_right - p.domain.y_left) / (p.sizeY - 1)),
                             2) / 5.0;

    if (solver.value("custom_integration_time_step", false)) {
        p.integrationTimeStep = solver.value("integration_time_step", computed_dt);
    } else {
        p.integrationTimeStep = computed_dt;
    }

    p.init_cond_from_file = solver.value("init_cond_from_file", false);
    p.init_cond_file_path = solver.value("init_cond_file_path", "");

    std::string ic_str = solver.value("initial_condition", "hyperbolic_tangent");

    if (ic_str == "hyperbolic_tangent")      p.init_condition = InitialCondition::HyperbolicTangent;
    else if (ic_str == "linear_by_parts")    p.init_condition = InitialCondition::LinearByParts;
    else if (ic_str == "constant_circle")    p.init_condition = InitialCondition::ConstantCircle;
    else if (ic_str == "constant_halves")    p.init_condition = InitialCondition::ConstantHalves;
    else if (ic_str == "stripe")             p.init_condition = InitialCondition::Stripe;
    else if (ic_str == "two_bumps")          p.init_condition = InitialCondition::TwoBumps;
    else if (ic_str == "star")               p.init_condition = InitialCondition::Star;
    else if (ic_str == "fourier_x")          p.init_condition = InitialCondition::FourierX;
    else if (ic_str == "fourier_y")          p.init_condition = InitialCondition::FourierY;
    else throw std::runtime_error("Unknown initial_condition in config: " + ic_str);

    auto problem = j.at("problem");
    p.alpha = problem.value("alpha", 1.0);
    p.beta  = problem.value("beta", 1.0);
    p.par_a = problem.value("a", 1.0);
    p.par_b = problem.value("b", 0.1);
    p.par_d = problem.value("D", 5e15);
    p.T = problem.value("T", 1200);
    p.ksi   = problem.value("ksi", 0.01);

    return p;
}

void Parameters::save_human_readable(const std::filesystem::path& filename) const {
    std::ofstream file(filename);
    if (!file) {
        throw std::runtime_error("Unable to open the file " + filename.string());
    }

    file << std::left << std::setw(24) << "Solver type:" << std::right << std::setw(28) << type << std::endl;
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
    file << std::left << std::setw(24) << "Alpha:" << std::right << std::setw(28) << alpha << std::endl;
    file << std::left << std::setw(24) << "Beta:"  << std::right << std::setw(28) << beta << std::endl;
    file << std::left << std::setw(24) << "Par_a:" << std::right << std::setw(28) << par_a << std::endl;
    file << std::left << std::setw(24) << "Par_b:" << std::right << std::setw(28) << par_b << std::endl;
    file << std::left << std::setw(24) << "Par_d:" << std::right << std::setw(28) << par_d << std::endl;
    file << std::left << std::setw(24) << "T:" << std::right << std::setw(28) << T << std::endl;
    file << std::left << std::setw(24) << "Ksi:"   << std::right << std::setw(28) << ksi << std::endl;
    file << std::left << std::setw(24) << "Model:" << std::right << std::setw(28) << static_cast<int>(model) << std::endl;
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