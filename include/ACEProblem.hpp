#pragma once

#include<cmath>
#include<iostream>
#include<fstream>
#include<vector>
#include<sstream>
#include<iomanip>
#include <filesystem>
#include <unistd.h>
#include <cassert>

#include "ODEProblem.hpp"

struct Domain
{
    double x_left;
    double x_right;
    double y_left;
    double y_right;
};

enum class MODEL
{
    MODEL_1 = 1,
    MODEL_2 = 2,
    MODEL_3 = 3
};

struct Parameters{
    double initial_time;
    double final_time;
    Domain domain;
    int sizeX;
    int sizeY;
    double timeStep;
    double integrationTimeStep;
    double alpha;
    double beta;
    double par_a;
    double ksi;
    MODEL model;
};


class ACEProblem : public ODEProblem
{
   public:

    ACEProblem(int sizeX,
               int sizeY,
               Domain domain,
               double alpha,
               double beta,
               double par_a,
               double ksi,
               MODEL model,
               std::string output_folder);
      
    int getDegreesOfFreedom();
      
    void getRightHandSide( const double& t, double* _u, double* fu );
      
    bool writeSolution( const double& t, int step, const double* u );

    void setInitialCondition( double* u );

    void set_dirichlet_boundary(double* _u, double* fu);

    double right_hand_side_at(double* _u, int i, int j);
    double laplace(double *u, int i, int j);
    double grad_norm(double *_u, int i, int j);

    double f_0(double *_u, int i, int j);

    double F(double *_u, int i, int j);

    protected:

    const int sizeX;
    const int sizeY;
    Domain domain;
    double hx;
    double hy;

    const double alpha;
    const double par_a;
    const double beta;
    const double ksi;

    const MODEL model;
    
    const std::string output_folder;
};