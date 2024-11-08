/*
 * File:   Euler.h
 * Author: oberhuber
 *
 * Created on February 25, 2016, 9:10 AM
 */

#pragma once

#include<algorithm>
#include "ODESolver.h"
#include "ODEProblem.h"

class Euler : public ODESolver
{
   public:

      Euler();

      bool setup( const int degreesOfFreedom );

      bool solve( const double integrationTimeStep,
                  const double stopTime,
                  double* time,
                  ODEProblem* problem,
                  double* u );

      ~Euler();

   protected:

      double* k;
};

