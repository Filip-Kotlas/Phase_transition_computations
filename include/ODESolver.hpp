/* 
 * File:   ODESolver.h
 * Author: oberhuber
 *
 * Created on February 25, 2016, 9:13 AM
 */

#pragma once

#include "ODEProblem.hpp"

class ODESolver
{
   public:
            
      virtual bool setup( const int degreesOfFreedom ) = 0;
      
      virtual bool solve( const double integrationTimeStep,
                          const double stopTime,
                          double* time,
                          ODEProblem* problem,
                          double* u ) = 0;
};

