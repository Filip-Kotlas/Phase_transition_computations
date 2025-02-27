/*
 * File:   ode-solve.h
 * Author: oberhuber
 *
 * Created on February 25, 2016, 8:48 AM
 */


#include <iostream>
#include <cmath>
#include <chrono>
#include <iomanip>
#include "ODEProblem.h"
#include "ODESolver.h"

bool solve( const double initialTime,
            const double finalTime,
            const double timeStep,
            const double integrationTimeStep,
            ODEProblem* problem,
            ODESolver* solver,
            double* u )
{
   solver->setup( problem->getDegreesOfFreedom() );
   const int timeStepsCount = std::ceil( std::max( 0.0, finalTime - initialTime ) / timeStep );
   double time( initialTime );
   problem->writeSolution( time, 0, u );
   auto start = std::chrono::high_resolution_clock::now();

   for( int k = 1; k <= timeStepsCount; k++ )
   {
      auto stepStart = std::chrono::high_resolution_clock::now();
      double currentTimeStep = std::min( timeStep, finalTime - time );
      if( !  solver->solve( integrationTimeStep,
                            time + currentTimeStep,  // stopTime
                            &time,
                            problem,
                            u ) )
         return false;

      problem->writeSolution( time, k, u );

      auto stepEnd = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> stepDuration = stepEnd - stepStart;
      
      double average_time_per_step = (std::chrono::duration<double>(stepEnd - start).count() / k);
      double remaining_time = average_time_per_step * (timeStepsCount - k);
      
      int hours = static_cast<int>(remaining_time) / 3600;
      int minutes = (static_cast<int>(remaining_time) % 3600) / 60;
      int seconds = remaining_time - (hours * 3600 + minutes * 60);
      
      std::cout << "Step " << k << " / " << timeStepsCount << " => " << std::fixed
                << std::setprecision(2) << ( double ) k / ( double ) timeStepsCount * 100.0 << "% ";
      std::cout << "     time remaining: " 
                << hours << "h " << minutes << "m " << seconds << "s"
                << std::endl;
   }
   std::cout << "Done." << std::endl;
   return true;
}

