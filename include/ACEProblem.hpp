#pragma once

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <iomanip>
#include <filesystem>
#include <unistd.h>
#include <cassert>
#include <algorithm>
#include <list>

#include "ODEProblem.hpp"
#include "constants.hpp"

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
    MODEL_3 = 3,
    MODEL_4 = 4
};

struct Parameters{
    std::string type;
    double initial_time;
    double final_time;
    Domain domain;
    int sizeX;
    int sizeY;
    int frame_num;
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
    * Phase and concentration access
    */
    double phase_at(const double* u, int i, int j){
        return u[j*sizeX + i];
    };
    double conc_at(const double* u, int i, int j){
        // Returns concentration in range (c_min, c_max)
        return u[sizeX * sizeY + j*sizeX + i];
    };

    /*
    * RHS computations
    */
    void getRightHandSide( const double& t, double* u, double* fu );
    double get_rhs_phase_at(const double* u, int i, int j);
    double get_rhs_concentration_at(const double &t, double* u, int i, int j);
    
    /*
    * Initial conditions
    */
    void setInitialCondition( double* u );
    void set_phase_initial_condition( double* u);
    void set_concentration_initial_condition( double* u);

    /*
    * Boundary conditions
    */
    void apply_boundary_condition(double* u, double* fu);
    void apply_phase_boundary_condition(double* u, double* fu);
    void apply_concentration_boundary_condition(double* u, double* fu);

    /*
    * Condition on concentration physical meaning
    */
    void apply_concentration_physical_condition(double* u);

    /*
    * Operators
    */
    double laplace(const double *u, int i, int j);
    double grad_norm(const double *_u, int i, int j);
    double div_D_grad_concentration(const double *u, int i, int j);
    double div_D_grad_phase(const double *u, int i, int j);
    double get_conc_diff_coef(const double *u, int i, int j);
    double get_phas_diff_coef(const double *u, int i, int j);

    double f_0(const double *_u, int i, int j);

    double grade_4_polynom(const double *u, int i, int j);
    double polynom_p(const double *u, int i, int j);
    double der_polynom_p(const double *u, int i, int j);

    double F(const double *u, int i, int j);
    double G(const double &t, double *u, int i, int j);

    double sec_deriv_of_g_w_resp_to_c(const double* u, int i, int j);
    double deriv_of_g_w_resp_to_c_and_p(const double* u, int i, int j);


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

    const double T = 1100;
    
    const MODEL model;
    
    const std::string output_folder;
};