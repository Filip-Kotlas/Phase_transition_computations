#include <iostream>
#include "Euler.h"

Euler::Euler()
{
   this-> k = 0;
}

bool Euler::setup( const int degreesOfFreedom )
{
   this->k = new double[ degreesOfFreedom ];
   if( ! this->k )
      return false;
   return true;
}

bool Euler::solve( const double integrationTimeStep,
                   const double stopTime,
                   double* time,
                   ODEProblem* problem,
                   double* u )
{
   const int dofs = problem->getDegreesOfFreedom();
   long int iteration( 0 );
   while( *time < stopTime )
   {
     const double tau = std::min( integrationTimeStep, stopTime - *time );
     problem->getRightHandSide( *time, u, this->k );
     for( int i = 0; i < dofs; i++ )
        u[ i ] +=  tau * k[ i ];
      *time += tau;
      iteration++;
      std::cout << "ITER: " << iteration << " \t tau = " << tau << " \t time= " << *time << "         \r " << std::flush;
   }
   std::cout << std::endl;
   return true;
}

Euler::~Euler()
{
   if( k ) delete[] k;
}
