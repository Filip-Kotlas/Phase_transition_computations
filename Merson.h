/*
 * File:   Merson.h
 * Author: oberhuber
 *
 * Created on February 26, 2016, 4:46 PM
 */

#pragma once

#include "ODESolver.h"
#include "ODEProblem.h"

class Merson : public ODESolver
{
   public:

      Merson();

      bool setup( int degreesOfFreedom );

      void setAdaptivity( const double& adaptivity );

      bool solve( const double integrationTimeStep,
                  const double stopTime,
                  double* time,
                  ODEProblem* problem,
                  double* u );

      ~Merson();

   protected:

      double *k1, *k2, *k3, *k4, *k5, *aux;

      double adaptivity;

};

