#pragma once
#include "types.hpp"
#include "constants.hpp"
#include "Parameters.hpp"

class InitialCondition {
public:

    InitialCondition() = default;

    InitialCondition(Parameters param)
    : type(param.init_condition), domain(param.domain), sizeX(param.sizeX), sizeY(param.sizeY), ksi(param.ksi), param(param) {
        r = (domain.x_right - domain.x_left)/6;
        r1 = r - 0.5*ksi;
        r2 = r1 + ksi;
        hx = (domain.x_right - domain.x_left) / (sizeX - 1);
        hy = (domain.y_right - domain.y_left) / (sizeY - 1);
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
    ICType type;
    Domain domain;
    Index sizeX;
    Index sizeY;
    Real ksi;
    Parameters param;
    

    Real r;
    Real r1;
    Real r2;
    Real hx;
    Real hy;
};