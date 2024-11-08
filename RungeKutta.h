/*
 * File:   Merson.h
 * Author: oberhuber
 *
 * Created on February 26, 2016, 4:46 PM
 */

#pragma once

#include "ODESolver.h"
#include "ODEProblem.h"

class RungeKutta : public ODESolver
{
   public:

      RungeKutta();

      bool setup( int degreesOfFreedom );

      bool solve( const double integrationTimeStep,
                  const double stopTime,
                  double* time,
                  ODEProblem* problem,
                  double* u );

      ~RungeKutta();

   protected:

      double *k1, *k2, *k3, *k4, *k5, *aux;

      double adaptivity;

};

