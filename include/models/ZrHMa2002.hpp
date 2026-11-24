#pragma once

// The hydride precipitation in the Zr-H alloy after Ma et al. (2002) (config: "model": "Zr-H_Ma2002"), with one
// orientation variant: the phase is the order parameter phi of the hydride (eta_1 in the article), the
// concentration is the hydrogen content C. The free energy (eq. 1 and 2 of the article, one variant)
//     F = int f(C, phi) + alpha/2 |grad phi|^2 + beta/2 |grad C|^2,
//     f = A1/2 (C - C1)^2 + A2/2 (C - C2) phi^2 - A3/4 phi^4 + A3/6 phi^6,
// and the equations (eq. 6: time-dependent Ginzburg-Landau for phi, Cahn-Hilliard for C)
//     d phi/dt = -L_p dF/dphi = L_p (alpha lap phi - A2 (C - C2) phi + A3 phi^3 - A3 phi^5),
//     d C/dt   = M lap mu,   mu = dF/dC = A1 (C - C1) + A2/2 phi^2 - beta lap C.
// Not included from the article: the elastic energy (eq. 5, which makes the hydrides needle-shaped), the
// noise terms (nucleation), the variants 2 and 3 with their cross terms A5-A7.
//
// The article has A4 phi^6/6 with its own A4, and its hydride is phi_article ~ 1.69 (sqrt(A3/A4)). Here
// phi = phi_article / k with k^2 = A3/A4 (of the article), so that the hydride is phi = 1 at C = C2 (and
// ~1 near it; f depends on phi^2 only, so -1 is a hydride too). The rescaling is exact (the same solution,
// divided by k) and makes the coefficients of phi^4/4 and phi^6/6 equal, so A4 is not a parameter here:
//     A2 = A2_article k^2, A3 = A3_article k^4 (= A4_article k^6), alpha = alpha_article k^2,
//     L_p = L_p_article / k^2; C1, C2, A1, beta, M unchanged.
// The article's values (units m k_B T and grid spacing 1) C1 = 0.05, C2 = 0.58, A1 = 60, A2 = -20, A3 = 10,
// A4 = 3.5, alpha = beta = 3.5, L_p = M = 0.4 give k^2 = 20/7 and A2 = -57.142857, A3 = 81.632653,
// alpha = 10, L_p = 0.14. The interface is only 2-3 grid spacings wide then: use h <= 0.5 (in those
// units); at h = 1 the bulk concentrations of a flat interface are 0.129 / 0.617 instead of the
// common tangent 0.104 / 0.58.
//
// Discretisation: 5-point Laplacians; the Cahn-Hilliard term is the 5-point Laplacian of mu, where mu is
// computed at the neighbours on the fly (a 13-point stencil for C). The equation needs two conditions on
// the boundary: dC/dn = 0 (the Neumann condition of the concentration, which check_parameters requires)
// and the no-flux condition dmu/dn = 0, in the same "copy the neighbour" form (mu on the boundary =
// mu at the neighbouring interior point). The sum of C over the interior points is then conserved.
// The explicit time step is limited by M (A1 k^2 + beta k^4), k^2 up to 4/hx^2 + 4/hy^2: about 0.01
// for the article's values and h = 1, and it falls like h^4 under refinement.
//
// The initial condition computed in code has the concentration c_init_hydride where the phase is 1 and
// c_init_matrix where it is 0 (the article: 0.58 and 0.1). The model has no interface width ksi and no
// anisotropy, so it cannot use the initial conditions that need them (HyperbolicTangent, WulffShape):
// start in main.cu refuses them.
//
// Its config section (all keys required, all floating-point numbers):
//     "Zr-H_Ma2002": { "C1": ..., "C2": ..., "A1": ..., "A2": ..., "A3": ...,
//                      "alpha": ..., "beta": ..., "L_p": ..., "M": ...,
//                      "c_init_hydride": ..., "c_init_matrix": ... }

#include <vector>
#include "external/json.hpp"

#include "models/Model.hpp"

class ZrHMa2002 : public Model
{
   public:

    // The name of the model in the config (the entry model).
    static constexpr const char* name = "Zr-H_Ma2002";

    // The physical parameters of this model (the config section parameters). Plain data,
    // copied to the GPU with the model.
    struct Params {
        Real C1;     // equilibrium hydrogen content of the matrix
        Real C2;     // equilibrium hydrogen content of the hydride
        Real A1;
        Real A2;
        Real A3;     // of phi^4/4 and phi^6/6 (the field is rescaled, see above)
        Real alpha;  // gradient energy coefficient of the phase (alpha_p in the article)
        Real beta;   // gradient energy coefficient of the concentration
        Real L_p;    // kinetic coefficient of the phase
        Real M;      // mobility of the hydrogen
        Real c_init_hydride;  // concentration of the initial condition where the phase is 1
        Real c_init_matrix;   // ... where it is 0

        // Reads the section strictly (ConfigSection.hpp); throws std::runtime_error on a wrong config.
        static Params load(const nlohmann::json& section);
        // The parameters for parameters.txt and parameters.tex.
        std::vector<ParameterEntry> entries() const;
        // The concentrations of the initial condition: c_init_hydride and c_init_matrix.
        InitialConcentrations initial_concentrations() const;
    };

    // The Cahn-Hilliard equation needs Neumann boundary conditions for the concentration; throws
    // std::runtime_error otherwise (called by start before the run folder is created).
    static void check_parameters(const Parameters& parameters);

    ZrHMa2002(const Parameters& parameters, const Params& own) : Model(parameters), own(own) {}

    __cuda_callable__
    void set_rhs_at(const VectorView& u, VectorView& fu, Index i, Index j) const;

   private:

    __cuda_callable__
    Real get_rhs_phase_at(const VectorView& u, Index i, Index j) const;
    __cuda_callable__
    Real get_rhs_concentration_at(const VectorView& u, Index i, Index j) const;

    // The chemical potential mu = dF/dC at (i, j); a point on the boundary is moved to the neighbouring
    // interior point (the no-flux condition, see above).
    __cuda_callable__
    Real chemical_potential_at(const VectorView& u, Index i, Index j) const;

    // 5-point Laplacians at an interior point.
    __cuda_callable__
    Real laplace_phase_at(const VectorView& u, Index i, Index j) const;
    __cuda_callable__
    Real laplace_conc_at(const VectorView& u, Index i, Index j) const;

    const Params own;
};
