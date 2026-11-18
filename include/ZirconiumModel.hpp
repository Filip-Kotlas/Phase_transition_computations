#pragma once

// The phase transition in the Zr-Nb alloy (force_term.type "zirconium"): the anisotropic
// Allen-Cahn equation for the phase, with the force term F = G_m_alpha - G_m_beta (the difference
// of the Gibbs energies, constants.hpp), coupled with the diffusion of the Nb concentration.
// See Model.hpp for what a model provides.

#include "Model.hpp"

class ZirconiumModel : public Model
{
   public:

    explicit ZirconiumModel(const Parameters& parameters) : Model(parameters) {}

    __cuda_callable__
    void set_rhs_at(const VectorView& u, VectorView& fu, Index i, Index j) const;
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
};
