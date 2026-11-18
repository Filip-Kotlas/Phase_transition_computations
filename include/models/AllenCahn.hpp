#pragma once

// The anisotropic Allen-Cahn equation for the phase with a prescribed force term F
// (model.type "allen_cahn", force_term "constant" or "radial"): the verification cases with
// analytic solutions (a shrinking circle or Wulff shape). The concentration is carried along but
// does not evolve. See Model.hpp for what a model provides.
//
// Its config section (all keys required):
//     "allen_cahn": { "alpha": 1, "a": 1, "b": 1, "ksi": 0.05, "A": 0.06, "m": 4, "theta_0": 0.0,
//                     "force_term": "constant", "force_size": 1.0 }   (or "radial", without force_size)

#include <vector>
#include "external/json.hpp"

#include "models/Model.hpp"

class AllenCahn : public Model
{
   public:

    // The name of the model in the config: model.type and the section model.allen_cahn.
    static constexpr const char* name = "allen_cahn";

    // The physical parameters of this model (the config section model.allen_cahn). Plain data,
    // copied to the GPU with the model.
    struct Params {
        // the anisotropic Allen-Cahn equation for the phase (ksi and anisotropy also go to the initial condition)
        Real alpha;
        Real par_a;   // a
        Real par_b;   // b
        Real ksi;
        Anisotropy anisotropy;   // A, m, theta_0
        FTType force_term;   // "constant": F = force_size, "radial": F = -2/r
        Real force_size;     // only for "constant"

        // Reads the section strictly (ConfigSection.hpp); throws std::runtime_error on a wrong config.
        static Params load(const nlohmann::json& section);
        // The parameters for parameters.txt and parameters.tex.
        std::vector<ParameterEntry> entries() const;
    };

    AllenCahn(const Parameters& parameters, const Params& own) : Model(parameters), own(own) {}

    __cuda_callable__
    void set_rhs_at(const VectorView& u, VectorView& fu, Index i, Index j) const;

    /*
    * The anisotropic Allen-Cahn equation for the phase (the same in both models)
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

    // The prescribed force term at (i, j).
    __cuda_callable__
    Real F(Index i, Index j) const;

   private:

    const Params own;
};
