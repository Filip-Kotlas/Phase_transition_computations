#pragma once

#include<cmath>
#include<iostream>
#include<fstream>
#include<vector>
#include<sstream>
#include<iomanip>
#include <filesystem>
#include <unistd.h>

#include "ODEProblem.h"


class ACEProblem : public ODEProblem
{
   public:

    ACEProblem( int sizeX, int sizeY, double alpha, double sigma, double ksi);
      
    int getDegreesOfFreedom();
      
    void getRightHandSide( const double& t, double* _u, double* fu );
      
    bool writeSolution( const double& t, int step, const double* u );

    void setInitialCondition( double* u );


    protected:

    const int sizeX;
    const int sizeY;
    double hx;
    double hy;

    double alpha;
    double sigma;
    double ksi;
};