#pragma once

// The discretised model: the right-hand side of the phase and concentration equations at the
// interior points and the boundary conditions, for the state vector u (layout in Grid.hpp).
// A Model is captured by value in the GPU lambdas of main.cu, so it holds only plain data
// (Grid, ModelParams) and every method used there is __cuda_callable__. It does not change after
// the construction, so all methods are const.

#include <cmath>
#include <algorithm>

#include "constants.hpp"
#include "Parameters.hpp"
#include "ModelParams.hpp"
#include "Grid.hpp"
#include "types.hpp"
#include "enums.hpp"

class Model
{
   public:

    explicit Model(const Parameters& parameters);

    Index getDegreesOfFreedom() const;

    /*
    * Phase and concentration access
    */
    __cuda_callable__
    Real phase_at(const VectorView& u, Index i, Index j) const {
        return u.getElement(grid.index(i, j));
    };
    __cuda_callable__
    Real conc_at(const VectorView& u, Index i, Index j) const {
        return u.getElement(grid.num_points() + grid.index(i, j));
    };

    /*
    * RHS computations
    */
    __cuda_callable__
    void set_rhs_at(const VectorView& u, VectorView& fu, Index i, Index j) const;
    __cuda_callable__
    Real get_rhs_phase_at(const VectorView& u, Index i, Index j) const;
    __cuda_callable__
    Real get_rhs_concentration_at(const VectorView& u, Index i, Index j) const;

    /*
    * Boundary conditions
    */
    __cuda_callable__
    void apply_boundary_condition_xdir(Index ind, VectorView fu) const;
    __cuda_callable__
    void apply_boundary_condition_ydir(Index ind, VectorView fu) const;
    __cuda_callable__
    void apply_phase_boundary_condition_xdir(Index ind, VectorView fu) const;
    __cuda_callable__
    void apply_phase_boundary_condition_ydir(Index ind, VectorView fu) const;
    __cuda_callable__
    void apply_concentration_boundary_condition_xdir(Index ind, VectorView fu) const;
    __cuda_callable__
    void apply_concentration_boundary_condition_ydir(Index ind, VectorView fu) const;

    /*
    * Operators
    */
    __cuda_callable__
    Real grad_p_1_forward(const VectorView& u, Index i, Index j) const;
    __cuda_callable__
    Real grad_p_2_forward(const VectorView& u, Index i, Index j) const;
    __cuda_callable__
    Real div_T0(const VectorView& u, Index i, Index j) const;
    __cuda_callable__
    Real T0_1(const Real grad_p_1, const Real grad_p_2) const;
    __cuda_callable__
    Real T0_2(const Real grad_p_1, const Real grad_p_2) const;
    __cuda_callable__
    Real div_D_grad_concentration(const VectorView& u, Index i, Index j) const;
    __cuda_callable__
    Real div_D_grad_phase(const VectorView& u, Index i, Index j) const;
    __cuda_callable__
    Real get_conc_diff_coef(const VectorView& u, Index i, Index j) const;
    __cuda_callable__
    Real get_phas_diff_coef(const VectorView& u, Index i, Index j) const;

    __cuda_callable__
    Real f_0(const VectorView& u, Index i, Index j) const;

    __cuda_callable__
    Real grade_4_polynom(const VectorView& u, Index i, Index j) const;
    __cuda_callable__
    Real polynom_p(const VectorView& u, Index i, Index j) const;
    __cuda_callable__
    Real der_polynom_p(const VectorView& u, Index i, Index j) const;

    __cuda_callable__
    Real F(const VectorView& u, Index i, Index j) const;

    __cuda_callable__
    Real sec_deriv_of_g_w_resp_to_c(const VectorView& u, Index i, Index j) const;
    __cuda_callable__
    Real deriv_of_g_w_resp_to_c_and_p(const VectorView& u, Index i, Index j) const;


    protected:

    const Grid grid;
    const ModelParams params;
};
