/* 
 * File:   ODEProblem.h
 * Author: oberhuber
 *
 * Created on February 5, 2018, 4:55 PM
 */
#pragma once

class ODEProblem
{
   public:
      
      virtual int getDegreesOfFreedom() = 0;
      
      virtual void getRightHandSide( const double& t, double* _u, double* fu ) = 0;
      
      virtual bool writeSolution( const double& t, int step, const double* u ) = 0;
};