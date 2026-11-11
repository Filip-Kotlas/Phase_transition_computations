#pragma once
#include <filesystem>
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

    // Helper functions: a simple hash and a uniform [0,1)
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

// Sets u (of length 2*sizeX*sizeY) to the phase and concentration computed by init_cond.
void set_init_cond_manually(Vector& u, const InitialCondition& init_cond, Index sizeX, Index sizeY);

// Sets u (of length 2*sizeX*sizeY) from a frame Result-*.npy of a sizeX x sizeY grid (see read_frame).
// Without scaling the frame is copied as it is. With scaling every second point of the frame
// (a grid with half the resolution) is placed into the centre of the domain, the rest is the beta phase.
// Prints the reason to std::cerr and returns false on failure.
bool set_init_cond_from_file(Vector& u, const std::filesystem::path& filename, Index sizeX, Index sizeY, bool scaling);