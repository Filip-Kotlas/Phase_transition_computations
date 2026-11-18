#pragma once

// The base class of the discretised models: what all of them share. Every model is a child class
// (AllenCahnModel.hpp, ZirconiumModel.hpp) that provides
//     Index getDegreesOfFreedom() const;                                   (inherited from here)
//     void set_rhs_at(const VectorView& u, VectorView& fu, Index i, Index j) const;   (interior points)
//     void apply_boundary_condition_xdir(Index ind, VectorView fu) const;  (inherited from here)
//     void apply_boundary_condition_ydir(Index ind, VectorView fu) const;  (inherited from here)
// main.cu calls them through the template run<ModelT>, so the model is chosen at compile time, not
// through virtual methods: a model is captured by value in the GPU lambdas, and an object with
// virtual methods copied from the host would carry a pointer to the host vtable onto the GPU.
// So a model holds only plain data (Grid, ModelParams) and every method used there is
// __cuda_callable__. It does not change after the construction, so all methods are const.
//
// All models so far have the same state vector u: the phase, then the concentration (layout in
// Grid.hpp), and the same boundary conditions.

#include "Parameters.hpp"
#include "ModelParams.hpp"
#include "Grid.hpp"
#include "types.hpp"
#include "enums.hpp"

class Model
{
   public:

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
    * The anisotropic Allen-Cahn equation for the phase
    */
    // The RHS of the phase equation at (i, j) for the value F of the force term there
    // (each model computes its own F).
    __cuda_callable__
    Real get_rhs_phase_at(const VectorView& u, Index i, Index j, Real F) const;

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
    Real f_0(const VectorView& u, Index i, Index j) const;
    __cuda_callable__
    Real grade_4_polynom(const VectorView& u, Index i, Index j) const;


    protected:

    // Only the child classes are constructed.
    explicit Model(const Parameters& parameters);

    const Grid grid;
    const ModelParams params;
};
