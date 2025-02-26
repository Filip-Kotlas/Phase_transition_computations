#include <cstdlib>
#include "ode-solve.h"
#include "RungeKutta.h"
#include "ODEProblem.h"
#include "ODESolver.h"
#include "ACEProblem.h"
#include <cassert>

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