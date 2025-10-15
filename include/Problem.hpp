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
#include <string>

#include "constants.hpp"
#include "Parameters.hpp"
#include "types.hpp"

class Problem
{
   public:

    Problem(Parameters param);
    
    Index getDegreesOfFreedom();
    bool writeSolution( const Real& t, Index step, const VectorView& u, const std::filesystem::path& output_folder);

    __cuda_callable__
    Index get_size_x() { return sizeX; };
    __cuda_callable__
    Index get_size_y() { return sizeY; };

    /*
    * Phase and concentration access
    */
    __cuda_callable__
    Real phase_at(const VectorView& u, Index i, Index j){
        return u.getElement(j*sizeX + i);
    };
    __cuda_callable__
    Real conc_at(const VectorView& u, Index i, Index j){
        return u.getElement(sizeX * sizeY + j*sizeX + i);
    };

    /*
    * RHS computations
    */
    __cuda_callable__
    void set_rhs_at(const VectorView& u, VectorView& fu, Index i, Index j);
    __cuda_callable__
    Real get_rhs_phase_at(const VectorView& u, Index i, Index j);
    __cuda_callable__
    Real get_rhs_concentration_at(const VectorView& u, Index i, Index j);
    
    /*
    * Initial conditions
    */
    void set_init_cond_manually( Vector &u, InitialCondition ic);
    void set_phase_initial_condition( Vector &u, InitialCondition ic);
    void set_concentration_initial_condition( Vector& u, InitialCondition ic);
    bool set_init_cond_from_file(Vector& u, const std::filesystem::path& filename);

    /*
    * Boundary conditions
    */
    __cuda_callable__
    void apply_boundary_condition_xdir(Index ind, VectorView u, VectorView fu);
    __cuda_callable__
    void apply_boundary_condition_ydir(Index ind, VectorView u, VectorView fu);
    __cuda_callable__
    void apply_phase_boundary_condition_xdir(Index ind, VectorView u, VectorView fu);
    __cuda_callable__
    void apply_phase_boundary_condition_ydir(Index ind, VectorView u, VectorView fu);
    __cuda_callable__
    void apply_concentration_boundary_condition_xdir(Index ind, VectorView u, VectorView fu);
    __cuda_callable__
    void apply_concentration_boundary_condition_ydir(Index ind, VectorView u, VectorView fu);

    // CPU
    /*
    void apply_boundary_condition(double *u, double *fu);
    void apply_phase_boundary_condition(double* u, double* fu);
    void apply_concentration_boundary_condition(double* u, double* fu);
    */

    /*
     * Condition on concentration physical meaning
     */
    //void apply_concentration_physical_condition(double *u);

    /*
    * Operators
    */
    __cuda_callable__
    Real laplace(const VectorView& u, Index i, Index j);
    __cuda_callable__
    Real div_D_grad_concentration(const VectorView& u, Index i, Index j);
    __cuda_callable__
    Real div_D_grad_phase(const VectorView& u, Index i, Index j);
    __cuda_callable__
    Real get_conc_diff_coef(const VectorView& u, Index i, Index j);
    __cuda_callable__
    Real get_phas_diff_coef(const VectorView& u, Index i, Index j);

    __cuda_callable__
    Real f_0(const VectorView& u, Index i, Index j);

    __cuda_callable__
    Real grade_4_polynom(const VectorView& u, Index i, Index j);
    __cuda_callable__
    Real polynom_p(const VectorView& u, Index i, Index j);
    __cuda_callable__
    Real der_polynom_p(const VectorView& u, Index i, Index j);

    __cuda_callable__
    Real F(const VectorView& u, Index i, Index j);
    //double G(const double &t, double *u, int i, int j);

    __cuda_callable__
    Real sec_deriv_of_g_w_resp_to_c(const VectorView& u, Index i, Index j);
    __cuda_callable__
    Real deriv_of_g_w_resp_to_c_and_p(const VectorView& u, Index i, Index j);


    protected:

    const Index sizeX;
    const Index sizeY;
    Domain domain;
    Real hx;
    Real hy;

    const Real alpha;
    const Real par_a;
    const Real par_b;
    const Real par_d;
    const Real ksi;

    const Real T;
    
    const MODEL model;
};