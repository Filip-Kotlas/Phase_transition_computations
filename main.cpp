#include <cstdlib>
#include "ode-solve.h"
#include "RungeKutta.h"
#include "ODEProblem.h"
#include "ODESolver.h"
#include "ACEProblem.h"
#include <windows.h>
#include <cassert>
#include <iostream>

bool directoryExists(const std::string& path) {
    DWORD attribs = GetFileAttributesA(path.c_str());
    return (attribs != INVALID_FILE_ATTRIBUTES && (attribs & FILE_ATTRIBUTE_DIRECTORY));
}

bool createDirectory(const std::string& path) {
    return CreateDirectoryA(path.c_str(), NULL) || GetLastError() == ERROR_ALREADY_EXISTS;
}

int main(int argc, char** argv)
{
    const double initialTime( 0.0 );
    const double finalTime( 0.30 );
    const Domain domain = {-1, 1, -1, 1};
    const int sizeX = 200;
    const int sizeY = 200;
    const double timeStep( 0.001 );
    const double integrationTimeStep(pow(std::min(  (domain.x_right - domain.x_left)/(sizeX-1),
                                                    (domain.y_right - domain.y_left)/(sizeY-1)),
                                         2)/4);
    const double alpha( 1 );
    const double beta( 1 );
    const double par_a( 10 );
    const double ksi( 0.01 );
    const MODEL model(MODEL::MODEL_1);


    std::string result_path = "Results";
    std::string info_path = result_path + "\\info";
    std::string setup_path = info_path + "\\setup.txt";

    if (!directoryExists(result_path)) {
        if (createDirectory(result_path)) {
            std::cout << "Created folder: " << result_path << std::endl;
        } else {
            std::cerr << "Error: Could not create folder!" << std::endl;
        }
    }
    if (!directoryExists(info_path)) {
        if (createDirectory(info_path)) {
            std::cout << "Created folder: " << info_path << std::endl;
        } else {
            std::cerr << "Error: Could not create folder!" << std::endl;
        }
    }

    
    std::fstream file;
    file.open(setup_path, std::fstream::out | std::fstream::trunc);
    if(!file)
    {
       std::cout << "Unable to open the file " << setup_path << std::endl;
       return false;
    }
    file << std::left << std::setw(24) << "Initial time:"  << std::right << std::setw(28) << initialTime << std::endl;
    file << std::left << std::setw(24) << "Final time:"    << std::right << std::setw(28) << finalTime << std::endl;
    std::ostringstream oss;
    oss << "[" 
        << "(" << std::fixed << std::setprecision(2) << domain.x_left << ", " 
        << domain.x_right << ")("
        << domain.y_left << ", " 
        << domain.y_right << ")]";
    std::string domain_str = oss.str();
    file << std::left << std::setw(24) << "Domain: "       << std::right << std::setw(28) << domain_str << std::endl;
    file << std::left << std::setw(24) << "SizeX:"         << std::right << std::setw(28) << sizeX << std::endl;
    file << std::left << std::setw(24) << "SizeY:"         << std::right << std::setw(28) << sizeY << std::endl;
    file << std::left << std::setw(24) << "Time step:"     << std::right << std::setw(28) << timeStep << std::endl;
    file << std::left << std::setw(24) << "Integration time step:" 
                            << std::right << std::setw(28) << integrationTimeStep << std::endl;
    file << std::left << std::setw(24) << "Alpha:"         << std::right << std::setw(28) << alpha << std::endl;
    file << std::left << std::setw(24) << "Beta:"          << std::right << std::setw(28) << beta << std::endl;
    file << std::left << std::setw(24) << "Par_a:"         << std::right << std::setw(28) << par_a << std::endl;
    file << std::left << std::setw(24) << "Ksi:"           << std::right << std::setw(28) << ksi << std::endl;
    file << std::left << std::setw(24) << "Model:"         << std::right << std::setw(28) << int(model);
    

    //Is not working. God knows why.
    /*
    std::filesystem::path result_path = "Results";
    std::filesystem::path info_path = result_path / "info";
    if(!std::filesystem::exists(info_path) || !std::filesystem::is_directory(info_path))
    {
        if(std::filesystem::create_directories(info_path))
            std::cout << "Created folder: " << info_path << std::endl;
        else
            std::cout << "Error: Could not create folder: " << info_path << std::endl;
    }
    */

    ACEProblem problem = ACEProblem(sizeX, sizeY, domain, alpha, beta, par_a, ksi, model);
    RungeKutta integrator;
    //integrator.setAdaptivity( adaptivity );

    double* u = new double[problem.getDegreesOfFreedom()];
    problem.setInitialCondition( u );
    problem.writeSolution( 0.0, 0, u );
   
    if( ! solve( initialTime,
            finalTime,
            timeStep,
            integrationTimeStep,
            &problem,
            &integrator,
            u ) )
    {
        delete[] u;
        return EXIT_FAILURE;
    }
    delete[] u;
    return EXIT_SUCCESS;
    
    return 0;
}