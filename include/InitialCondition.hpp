#pragma once
#include "types.hpp"
#include "constants.hpp"
#include "Parameters.hpp"
#include "enums.hpp"

class InitialCondition {
public:

    InitialCondition() = default;

    InitialCondition(Parameters param)
    : param(param) {
        init_cond_radius = param.init_cond_radius;
        slope_width = param.init_cond_slope_width;
        hx = (param.domain.x_right - param.domain.x_left) / (param.sizeX - 1);
        hy = (param.domain.y_right - param.domain.y_left) / (param.sizeY - 1);
    }

    __cuda_callable__
    Real get_phase(Index ind) const;

    __cuda_callable__
    Real get_concentration(Index index) const;

    // Pomocné funkce: jednoduchý hash a uniformní [0,1)
    __cuda_callable__ static
    uint32_t wang_hash(uint32_t x);

    __cuda_callable__ static
    Real rnd01(uint32_t idx, uint32_t seed);

    __cuda_callable__
    Real psi(const Real theta) const;

    __cuda_callable__
    Real der_psi(const Real theta) const;


private:
    Parameters param;
    

    Real init_cond_radius;
    Real slope_width;
    Real hx;
    Real hy;
};