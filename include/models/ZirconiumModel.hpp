#pragma once

// The phase transition in the Zr-Nb alloy (model.type "zirconium"): the anisotropic
// Allen-Cahn equation for the phase, with the force term F = G_m_alpha - G_m_beta (the difference
// of the Gibbs energies, constants.hpp), coupled with the diffusion of the Nb concentration.
// See Model.hpp for what a model provides.
//
// Its config section (all keys required):
//     "zirconium": { "alpha": 1, "a": 1, "b": 3.53553391, "ksi": 0.075, "A": 0.06, "m": 4,
//                    "theta_0": -0.392699, "d": 1.6e15, "T": 1000 }

#include <vector>
#include "external/json.hpp"

#include "models/Model.hpp"

class ZirconiumModel : public Model
{
   public:

    // The name of the model in the config: model.type and the section model.zirconium.
    static constexpr const char* name = "zirconium";

    // The physical parameters of this model (the config section model.zirconium). Plain data,
    // copied to the GPU with the model.
    struct Params {
        // the anisotropic Allen-Cahn equation for the phase (ksi and anisotropy also go to the initial condition)
        Real alpha;
        Real par_a;   // a
        Real par_b;   // b
        Real ksi;
        Anisotropy anisotropy;   // A, m, theta_0
        Real d;   // the scaling of the diffusion
        Real T;   // the temperature [K]

        // Reads the section strictly (ConfigSection.hpp); throws std::runtime_error on a wrong config.
        static Params load(const nlohmann::json& section);
        // The parameters for parameters.txt and parameters.tex.
        std::vector<ParameterEntry> entries() const;
    };

    ZirconiumModel(const Parameters& parameters, const Params& own) : Model(parameters), own(own) {}

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

    __cuda_callable__
    Real get_rhs_concentration_at(const VectorView& u, Index i, Index j) const;

    // The force term at (i, j).
    __cuda_callable__
    Real F(const VectorView& u, Index i, Index j) const;

    /*
    * Diffusion of the concentration
    */
    __cuda_callable__
    Real div_D_grad_concentration(const VectorView& u, Index i, Index j) const;
    __cuda_callable__
    Real div_D_grad_phase(const VectorView& u, Index i, Index j) const;
    __cuda_callable__
    Real get_conc_diff_coef(const VectorView& u, Index i, Index j) const;
    __cuda_callable__
    Real get_phas_diff_coef(const VectorView& u, Index i, Index j) const;

    __cuda_callable__
    Real polynom_p(const VectorView& u, Index i, Index j) const;
    __cuda_callable__
    Real der_polynom_p(const VectorView& u, Index i, Index j) const;

    __cuda_callable__
    Real sec_deriv_of_g_w_resp_to_c(const VectorView& u, Index i, Index j) const;
    __cuda_callable__
    Real deriv_of_g_w_resp_to_c_and_p(const VectorView& u, Index i, Index j) const;

   private:

    const Params own;
};
