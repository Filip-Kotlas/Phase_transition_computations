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

    ACEProblem( int sizeX, int sizeY, double alpha, double sigma, double ksi, double T);
      
    int getDegreesOfFreedom();
      
    void getRightHandSide( const double& t, double* _u, double* fu );
      
    bool writeSolution( const double& t, int step, const double* u );

    void setInitialCondition( double* u );

    double laplaceD(double *u, int i, int j);

    double f_0(double *_u, int i, int j);

    double F(double *_u, int i, int j);

    protected:

    const int sizeX;
    const int sizeY;
    double hx;
    double hy;

    const double alpha;
    const double sigma;
    const double ksi;
    const double T;

    //parameters
    static const double delta_0 = 5e-9;
    static const double b = 3.23e-10;
    
};