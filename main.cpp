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
    const double finalTime( 0.1 );
    const double timeStep( 0.0001 );
    const double integrationTimeStep( 0.00001 );
    const int sizeX( 101 );
    const int sizeY( 101 );
    const double alpha( 3 );
    const double sigma( 1 );
    const double ksi ( 0.25 );
    
    ACEProblem problem = ACEProblem(sizeX, sizeY, alpha, sigma, ksi);
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