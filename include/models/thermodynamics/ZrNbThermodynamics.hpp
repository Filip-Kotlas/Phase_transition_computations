#pragma once

// The thermodynamic data of the Zr-Nb alloy (CALPHAD-style Gibbs energies, diffusivities, mobilities),
// used by ZrNbChoud2013. The values that depend only on the temperature T are computed once, on the
// host in double, by the constructor (the model builds it from its T before the time loop). It is
// plain data: the model holds it, so it goes to the GPU with the model, where the functions of the
// concentration are evaluated.

#include "core/types.hpp"

struct ZrNbThermodynamics
{
    // The equilibrium concentrations of Nb in the two phases (used by the initial condition).
    static constexpr Real c_init_alpha = 0.007;
    static constexpr Real c_init_beta = 0.025;

    // The interaction parameter of the alpha phase [J/mol] (it does not depend on T).
    static constexpr Real L_0_alpha = 24411;

    // The values at the temperature T [K]:
    Real RT;   // R*T [J/mol]
    // the Gibbs energies of the pure elements in the two phases [J/mol]
    Real G_Nb_alpha_0, G_Zr_alpha_0, G_Nb_beta_0, G_Zr_beta_0;
    // the interaction parameters of the beta phase [J/mol]
    Real L_0_beta, L_0_i_beta;
    // the mobilities of Nb in the two phases
    Real M_Nb_alpha, M_Nb_beta;

    // Computes the values above (host only).
    explicit ZrNbThermodynamics(double T);

    // The molar Gibbs energies of the two phases at the Nb concentration c [J/mol].
    __cuda_callable__ Real G_m_alpha(Real c) const;
    __cuda_callable__ Real G_m_beta(Real c) const;
};
