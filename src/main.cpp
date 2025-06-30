#include <cstdlib>
#include <windows.h>
#include <cassert>
#include <iostream>
#include <vector>

#include <json.hpp>

#include "ode-solve.hpp"
#include "Merson.hpp"
#include "ODEProblem.hpp"
#include "ODESolver.hpp"
#include "ACEProblem.hpp"


bool directoryExists(const std::string& path) {
    DWORD attribs = GetFileAttributesA(path.c_str());
    return (attribs != INVALID_FILE_ATTRIBUTES && (attribs & FILE_ATTRIBUTE_DIRECTORY));
}

bool createDirectory(const std::string& path) {
    return CreateDirectoryA(path.c_str(), NULL) || GetLastError() == ERROR_ALREADY_EXISTS;
}

Parameters get_parameters() {
    Parameters par;
    std::ifstream file("config\\config.json");
    if (!file) {
        throw std::runtime_error("Nelze otevřít soubor config\\config.json");
    }
    nlohmann::json config;
    file >> config;

    try{
        if (!config.contains("solver") || !config.contains("problem")) {
            throw std::runtime_error("Chybí sekce 'solver' nebo 'problem' v JSON.");
        }
        nlohmann::json solver = config["solver"];
        par.initial_time = solver.value("initial_time", 0.0);
        par.final_time = solver.value("final_time", 0.30);
        par.sizeX = solver.value("sizeX", 200);
        par.sizeY = solver.value("sizeY", 200);
        par.timeStep = (par.final_time - par.initial_time) / 100;
        
        int model_value = solver.value("model", -1);
        if( !(model_value == 1 || model_value == 2 || model_value == 3 || model_value == 4)) {
            throw std::runtime_error("Neplatná hodnota 'model' v JSON.");
        }
        par.model = static_cast<MODEL>(model_value);
        
        if (!solver.contains("domain")) {
            throw std::runtime_error("Chybí 'domain' v sekci 'solver'.");
        }
        
        nlohmann::json domain = solver["domain"];
        par.domain = {
            domain.value("x_left", -1.0),
            domain.value("x_right", 1.0),
            domain.value("y_left", -1.0),
            domain.value("y_right", 1.0)
        };

        par.integrationTimeStep = pow(std::min((par.domain.x_right - par.domain.x_left)/(par.sizeX-1),
                                               (par.domain.y_right - par.domain.y_left)/(par.sizeY-1)),
                                      2)/4;

        nlohmann::json problem = config["problem"];
        par.alpha = problem.value("alpha", 1.0);
        par.beta = problem.value("beta", 1.0);
        par.par_a = problem.value("a", 1.0);
        par.ksi = problem.value("ksi", 0.01);
    }
    catch(const std::exception& e) {
        std::cout << "Chyba při načítání JSON: " << e.what() << std::endl;
        throw;
    }

    return par;
}

void create_folders(std::vector<std::string> folder_paths) {
    for(const std::string& folder : folder_paths)
    {
        if (!directoryExists(folder)) {
            if (createDirectory(folder)) {
                std::cout << "Created folder: " << folder << std::endl;
            } else {
                std::cout << "Error: Could not create folder" << folder << "!" << std::endl;
            }
        }
    }
}

void save_parameters(std::string file_path, Parameters param) {
    std::fstream file;
    file.open(file_path, std::fstream::out | std::fstream::trunc);
    if(!file)
    {
       std::cout << "Unable to open the file " << file_path << std::endl;
       return;
    }
    file << std::left << std::setw(24) << "Initial time:"  << std::right << std::setw(28) << param.initial_time << std::endl;
    file << std::left << std::setw(24) << "Final time:"    << std::right << std::setw(28) << param.final_time << std::endl;
    std::ostringstream oss;
    oss << "[" 
        << "(" << std::fixed << std::setprecision(2) << param.domain.x_left << ", " 
        << param.domain.x_right << ")("
        << param.domain.y_left << ", " 
        << param.domain.y_right << ")]";
    std::string domain_str = oss.str();
    file << std::left << std::setw(24) << "Domain: "       << std::right << std::setw(28) << domain_str << std::endl;
    file << std::left << std::setw(24) << "SizeX:"         << std::right << std::setw(28) << param.sizeX << std::endl;
    file << std::left << std::setw(24) << "SizeY:"         << std::right << std::setw(28) << param.sizeY << std::endl;
    file << std::left << std::setw(24) << "Time step:"     << std::right << std::setw(28) << param.timeStep << std::endl;
    file << std::left << std::setw(24) << "Integration time step:" 
                            << std::right << std::setw(28) << param.integrationTimeStep << std::endl;
    file << std::left << std::setw(24) << "Alpha:"         << std::right << std::setw(28) << param.alpha << std::endl;
    file << std::left << std::setw(24) << "Beta:"          << std::right << std::setw(28) << param.beta << std::endl;
    file << std::left << std::setw(24) << "Par_a:"         << std::right << std::setw(28) << param.par_a << std::endl;
    file << std::left << std::setw(24) << "Ksi:"           << std::right << std::setw(28) << param.ksi << std::endl;
    file << std::left << std::setw(24) << "Model:"         << std::right << std::setw(28) << int(param.model);
}

void copy_config_file(std::string config_path, std::string output_path)
{
    std::fstream output_file;
    output_file.open(output_path, std::fstream::out | std::fstream::trunc);
    if(!output_file)
    {
       std::cout << "Unable to open the file " << output_path << std::endl;
       return;
    }

    std::fstream input_file;
    input_file.open(config_path, std::fstream::in);
    if(!input_file)
    {
       std::cout << "Unable to open the file " << config_path << std::endl;
       return;
    }

    output_file << input_file.rdbuf();
}

int main(int argc, char** argv)
{
    Parameters parameters = get_parameters();

    std::string result_path = "Results";
    if(argc == 2)
        result_path = argv[1];
    std::string info_path = result_path + "\\info";
    std::string setup_path = info_path + "\\parameters.txt";
    std::string calc_path = result_path + "\\calculations";
    std::string config_path = info_path + "\\config.json";

    create_folders({result_path, info_path, calc_path});
    save_parameters(setup_path, parameters);
    copy_config_file("config\\config.json", config_path);

    
    int arr[11] = {882, 911, 949, 993, 1100, 1149, 1350, 1450, 1545, 1645, 1758};
    for(int i = 0; i < 11; i++)
    {
        std::cout << arr[i] << ": " << constants::D_Nb_beta(arr[i] + 273.15) << ", " << constants::D_Nb_alpha(arr[i] + 273.15) << std::endl;
    }
    return 0;

    ACEProblem problem = ACEProblem(parameters.sizeX,
                                    parameters.sizeY,
                                    parameters.domain,
                                    parameters.alpha,
                                    parameters.beta,
                                    parameters.par_a,
                                    parameters.ksi,
                                    parameters.model,
                                    calc_path);
    Merson integrator;

    double* u = new double[problem.getDegreesOfFreedom()];
    problem.setInitialCondition( u );
    problem.writeSolution( 0.0, 0, u );
   
    if( ! solve(parameters.initial_time,
                parameters.final_time,
                parameters.timeStep,
                parameters.integrationTimeStep,
                &problem,
                &integrator,
                u))
    {
        delete[] u;
        return EXIT_FAILURE;
    }
    delete[] u;
    return EXIT_SUCCESS;
    
    return 0;
}