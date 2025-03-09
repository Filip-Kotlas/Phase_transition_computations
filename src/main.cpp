#include <cstdlib>
#include "ode-solve.hpp"
#include "RungeKutta.hpp"
#include "ODEProblem.hpp"
#include "ODESolver.hpp"
#include "ACEProblem.hpp"
#include <windows.h>
#include <cassert>
#include <iostream>
#include <vector>

bool directoryExists(const std::string& path) {
    DWORD attribs = GetFileAttributesA(path.c_str());
    return (attribs != INVALID_FILE_ATTRIBUTES && (attribs & FILE_ATTRIBUTE_DIRECTORY));
}

bool createDirectory(const std::string& path) {
    return CreateDirectoryA(path.c_str(), NULL) || GetLastError() == ERROR_ALREADY_EXISTS;
}

Parameters get_parameters() {
    Parameters par;
    par.initial_time = 0.0;
    par.final_time = 300;
    par.domain = {-1, 1, -1, 1};
    par.sizeX = 200;
    par.sizeY = 200;
    par.timeStep =  0.1;
    par.integrationTimeStep = pow(std::min((par.domain.x_right
                                                   -par.domain.x_left)
                                                   /(par.sizeX-1),
                                                  (par.domain.y_right
                                                   -par.domain.y_left)
                                                   /(par.sizeY-1)),
                                         2)/4;
    par.integrationTimeStep = 0.1;
    par.alpha = 1;
    par.beta = 1;
    par.par_a = 1;
    par.ksi = 0.01;
    par.model = MODEL::MODEL_3;

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

int main(int argc, char** argv)
{
    Parameters parameters = get_parameters();

    std::string result_path = "Results";
    if(argc == 2)
        result_path = argv[1];
    std::string info_path = result_path + "\\info";
    std::string setup_path = info_path + "\\setup.txt";
    std::string calc_path = result_path + "\\calculations";

    create_folders({result_path, info_path, calc_path});
    save_parameters(setup_path, parameters);

    ACEProblem problem = ACEProblem(parameters.sizeX,
                                    parameters.sizeY,
                                    parameters.domain,
                                    parameters.alpha,
                                    parameters.beta,
                                    parameters.par_a,
                                    parameters.ksi,
                                    parameters.model,
                                    calc_path);
    RungeKutta integrator;
    //integrator.setAdaptivity( adaptivity );

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