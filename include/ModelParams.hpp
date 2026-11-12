#pragma once

// The parameters of the model that the GPU code needs (Model, InitialCondition). Plain data only
// (numbers and enums), so that copying them to the GPU by value is well defined; Parameters, which
// also holds strings and the settings of the run, stays on the host. Built by Parameters::model_params().

#include <cmath>

#include "types.hpp"
#include "enums.hpp"

// The anisotropy psi(theta) = 1 + A cos(m (theta - theta_0)) of the interface energy.
struct Anisotropy {
    Real A;
    int m;
    Real theta_0;

    __cuda_callable__
    Real psi(const Real theta) const { return 1 + A*cos(m*(theta - theta_0)); }

    // d psi / d theta
    __cuda_callable__
    Real der_psi(const Real theta) const { return -A*m*sin(m*(theta - theta_0)); }
};

struct ModelParams {
    // phase and concentration equations (the config section "problem")
    Real alpha;
    Real par_a;
    Real par_b;
    Real par_d;
    Real T;
    Real ksi;
    // force term
    FTType force_term_type;
    Real force_term_size;
    Anisotropy anisotropy;
    // boundary conditions
    BCType bc_phase_x;
    BCType bc_phase_y;
    BCType bc_conc_x;
    BCType bc_conc_y;
};
