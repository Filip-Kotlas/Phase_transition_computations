#include <cstdlib>
#include "ode-solve.h"
#include "Merson.h"
#include "ODEProblem.h"
#include "ODESolver.h"
#include "ACEProblem.h"

int main(int argc, char** argv)
{
    const double initialTime( 0.0 );
    const double finalTime( 0.05 );
    const double timeStep( 0.0001 );
    const double integrationTimeStep( 0.0001 );
    const double adaptivity( 0 );
    const int sizeX( 100 );
    const int sizeY( 100 );
    const double alpha( 0.1 );
    const double sigma( 0.1 );
    const double ksi ( 0.01 );
    
    ACEProblem problem(sizeX, sizeY, alpha, sigma, ksi);
    Merson integrator;
    integrator.setAdaptivity( adaptivity );

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