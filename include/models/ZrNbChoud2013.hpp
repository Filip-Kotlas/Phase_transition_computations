#pragma once

// The phase transition in the Zr-Nb alloy (config: "model": "Zr-Nb_Choud2013"): the anisotropic
// Allen-Cahn equation for the phase, with the force term F = G_m_alpha - G_m_beta (the difference
// of the molar Gibbs energies of the two phases), coupled with the diffusion of the Nb concentration.
// The Gibbs energies come from a CALPHAD database in the TDB format (io/TdbDatabase.hpp), at the
// temperature T, with the concentration c = the mole fraction of Nb (BinaryPhase.hpp); the mobilities of
// Nb, which such databases do not hold, are in ZrNbChoud2013.cu. See Model.hpp for what a model provides.
//
// Its config section parameters (all keys required; database relative to the folder the program runs
// in, the repo root; the phases are the database's names of the alpha (phase 1) and beta (phase 0) phase):
//     "parameters": { "alpha": 1, "a": 1, "b": 3.53553391, "ksi": 0.075, "A": 0.06, "m": 4,
//                     "theta_0": -0.392699, "d": 1.6e15, "T": 1000,
//                     "database": "databases/NbZr-91Fer-LB.tdb", "alpha_phase": "HCP_A3", "beta_phase": "BCC_A2" }

#include <vector>
#include "external/json.hpp"

#include "models/Model.hpp"
#include "models/BinaryPhase.hpp"

class ZrNbChoud2013 : public Model
{
   public:

    // The name of the model in the config (the entry model).
    static constexpr const char* name = "Zr-Nb_Choud2013";

    // The physical parameters of this model (the config section parameters). Plain data,
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
        // The Gibbs energies of the two phases at T, from the database (the config entries database,
        // alpha_phase and beta_phase; Parameters::model_text_entries lists them in parameters.txt).
        // Nb is substitutional, so the concentration is the mole fraction of Nb. The phases have at most
        // the Redlich-Kister terms L_0 and L_1 (both Nb-Zr databases; the reader refuses a higher order).
        SubstitutionalPhase< 2 > alpha_phase, beta_phase;

        // Reads the section strictly (ConfigSection.hpp) and the database; throws std::runtime_error on a
        // wrong config or database.
        static Params load(const nlohmann::json& section);
        // The parameters for parameters.txt and parameters.tex (the numbers; the texts are added by
        // Parameters::model_text_entries).
        std::vector<ParameterEntry> entries() const;
        // The concentrations of the initial condition: the equilibrium ones of the two phases
        // (c_init_alpha where the phase is 1, c_init_beta where it is 0).
        InitialConcentrations initial_concentrations() const;
    };

    // The equilibrium concentrations of Nb in the two phases (used by the initial condition).
    static constexpr Real c_init_alpha = 0.007;
    static constexpr Real c_init_beta = 0.025;

    // Also computes the mobilities of Nb at the temperature T, once, on the host.
    ZrNbChoud2013(const Parameters& parameters, const Params& own);

    __cuda_callable__
    void set_rhs_at(const VectorView& u, VectorView& fu, Index i, Index j) const;

    /*
    * The anisotropic Allen-Cahn equation for the phase (the same in both models)
    */
    // The RHS of the phase equation at (i, j). The force term F (G_m_alpha - G_m_beta) is computed
    // inside, from the concentration there (AllenCahn's copy takes its prescribed F as an argument).
    __cuda_callable__
    Real get_rhs_phase_at(const VectorView& u, Index i, Index j) const;

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
    // The mobilities of Nb in the two phases at own.T (declared after own, so they are initialized after it).
    const Real M_Nb_alpha, M_Nb_beta;
};
