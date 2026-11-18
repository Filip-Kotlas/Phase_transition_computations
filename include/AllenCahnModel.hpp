#pragma once

// The anisotropic Allen-Cahn equation for the phase with a prescribed force term F
// (force_term.type "constant" or "radial"): the verification cases with analytic solutions
// (a shrinking circle or Wulff shape). The concentration is carried along but does not evolve.
// See Model.hpp for what a model provides.

#include "Model.hpp"

class AllenCahnModel : public Model
{
   public:

    explicit AllenCahnModel(const Parameters& parameters) : Model(parameters) {}

    __cuda_callable__
    void set_rhs_at(const VectorView& u, VectorView& fu, Index i, Index j) const;

    // The prescribed force term at (i, j).
    __cuda_callable__
    Real F(Index i, Index j) const;
};
