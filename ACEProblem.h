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

struct Domain
{
    double x_left;
    double x_right;
    double y_left;
    double y_right;
};

class ACEProblem : public ODEProblem
{
   public:

    ACEProblem( int sizeX, int sizeY, Domain domain, double alpha, double sigma, double ksi);
      
    int getDegreesOfFreedom();
      
    void getRightHandSide( const double& t, double* _u, double* fu );
      
    bool writeSolution( const double& t, int step, const double* u );

    void setInitialCondition( double* u );

    void set_dirichlet_boundary(double* _u, double* fu);

    double laplace(double *u, int i, int j);

    double f_0(double *_u, int i, int j);

    double F(double *_u, int i, int j);

    protected:

    const int sizeX;
    const int sizeY;
    Domain domain;
    double hx;
    double hy;

    const double alpha;
    const double sigma;
    const double ksi;    
};