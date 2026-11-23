#pragma once

// The base class of the discretised models: what all of them share. Every model is a child class
// (AllenCahn.hpp, ZrNbChoud2013.hpp, ZrHMa2002.hpp) that provides
//     static constexpr const char* name;   (model.type in the config, and its section model.<name>)
//     struct Params;                       (its physical parameters, plain data, with
//         static Params load(const nlohmann::json& section);   (reads model.<name>, see ConfigSection.hpp)
//         std::vector<ParameterEntry> entries() const;         (for parameters.txt and .tex)
//         InitialConcentrations initial_concentrations() const; (the concentrations of the initial
//             condition where the phase is 1 and 0, see Parameters.hpp; from constants or the config))
//     static void check_parameters(const Parameters& parameters);  (inherited from here, does nothing;
//         a model that cannot work with some run settings hides it with its own, which throws
//         std::runtime_error; start calls it before the run folder is created)
//     ModelT(const Parameters& parameters, const Params& own);
//     Index getDegreesOfFreedom() const;                                   (inherited from here)
//     void set_rhs_at(const VectorView& u, VectorView& fu, Index i, Index j) const;   (interior points)
//     void apply_boundary_condition_xdir(Index ind, VectorView fu) const;  (inherited from here)
//     void apply_boundary_condition_ydir(Index ind, VectorView fu) const;  (inherited from here)
// main.cu calls them through the templates start<ModelT> and run<ModelT>, so the model is chosen at
// compile time, not through virtual methods: a model is captured by value in the GPU lambdas, and an
// object with virtual methods copied from the host would carry a pointer to the host vtable onto the GPU.
// So a model holds only plain data (Grid, BoundaryConditions, its Params) and every method used
// there is __cuda_callable__. It does not change after the construction, so all methods are const.
//
// All models so far have the same state vector u: the phase, then the concentration (layout in
// Grid.hpp), and the same boundary conditions (BoundaryConditions, from the config section bound_cond).
// Everything else (the equations and their physical parameters, from the model's own config section)
// is in the child classes; the anisotropic Allen-Cahn operator for the phase is a copy in each of them
// that uses it (AllenCahn, ZrNbChoud2013).

#include "config/Parameters.hpp"
#include "core/Grid.hpp"
#include "core/types.hpp"
#include "config/enums.hpp"

class Model
{
   public:

    Index getDegreesOfFreedom() const;

    // Checks the run settings (e.g. the boundary conditions) that the model cannot work with. Nothing
    // here; a child class hides it with its own if it needs a check (see the list at the top).
    static void check_parameters(const Parameters&) {}

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

    protected:

    // Only the child classes are constructed.
    explicit Model(const Parameters& parameters);

    const Grid grid;
    const BoundaryConditions bound_cond;
};
