#include <cstdlib>
#include "ode-solve.h"
#include "RungeKutta.h"
#include "Euler.h"
#include "ODEProblem.h"
#include "ODESolver.h"
#include "ACEProblem.h"
#include <cassert>

int main(int argc, char** argv)
{
    const double initialTime( 0.0 );
    const double finalTime( 0.15 );
    const double timeStep( 0.001 );
    const double integrationTimeStep( 0.00002 );
    const Domain domain = {-1, 1, -1, 1};
    const int sizeX = 200;
    const int sizeY = 200;
    const double alpha( 1 );
    const double sigma( 1 );
    const double ksi ( 0.01 );
    
    ACEProblem problem = ACEProblem(sizeX, sizeY, domain, alpha, sigma, ksi);
    RungeKutta integrator;
    //integrator.setAdaptivity( adaptivity );

    double* u = new double[ problem.getDegreesOfFreedom() ];
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