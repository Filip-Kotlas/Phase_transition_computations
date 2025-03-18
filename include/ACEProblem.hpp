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
#include <algorithm>

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
    bool writeSolution( const double& t, int step, const double* u );

    /*
    * RHS computations
    */
    void getRightHandSide( const double& t, double* u, double* fu );
    double get_rhs_phase_at(double* u, int i, int j);
    double get_rhs_concentration_at(double* u, int i, int j);
    
    /*
    * Initial conditions
    */
    void setInitialCondition( double* u );
    void set_phase_initial_condition( double* u);
    void set_concentration_initial_condition( double* u);

    /*
    * Boundary conditions
    */
    void apply_boundary_condition(double* _u, double* fu);
    void apply_phase_boundary_condition(double* _u, double* fu);
    void apply_concentration_boundary_condition(double* _u, double* fu);

    double laplace(double *u, int i, int j);
    double grad_norm(double *_u, int i, int j);
    double div_D_grad_concentration(double *u, int i, int j);
    double get_diffusion_coef(double *u, int i, int j);

    double f_0(double *_u, int i, int j);

    double F(double *u, int i, int j);
    double G(double *u, int i, int j);

    double get_M_phi_tilde();
    double get_epsilon_phi_tilde();
    double get_G_alpha_tilde(const double *_u, int i, int j);
    double get_G_beta_tilde(const double *_u, int i, int j);
    double get_w_tilde();
    
    double get_p_prime(double *_u, int i, int j);
    double get_q_prime(double *_u, int i, int j);

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
    
    const int T = 1200;

    const MODEL model;
    
    const std::string output_folder;
};