#pragma once
#include <filesystem>
#include "types.hpp"
#include "constants.hpp"
#include "Parameters.hpp"
#include "ModelParams.hpp"
#include "Grid.hpp"
#include "enums.hpp"

// The initial condition computed in code (init_condition in the config). Captured by value in the
// GPU lambda of set_init_cond_manually, so it holds only plain data.
class InitialCondition {
public:

    explicit InitialCondition(const Parameters& parameters)
    :  init_condition(parameters.init_condition),
       init_cond_radius(parameters.init_cond_radius),
       slope_width(parameters.init_cond_slope_width),
       ksi(parameters.ksi),
       anisotropy(parameters.model_params().anisotropy),
       grid(parameters.grid())
    {}

    __cuda_callable__
    Real get_phase(Index ind) const;

    __cuda_callable__
    Real get_concentration(Index index) const;

    // Helper functions: a simple hash and a uniform [0,1)
    __cuda_callable__ static
    uint32_t wang_hash(uint32_t x);

    __cuda_callable__ static
    Real rnd01(uint32_t idx, uint32_t seed);

private:
    ICType init_condition;
    Real init_cond_radius;
    Real slope_width;
    Real ksi;                 // interface width, for HyperbolicTangent
    Anisotropy anisotropy;    // for WulffShape
    Grid grid;
};

// Sets u (of length 2*grid.num_points()) to the phase and concentration computed by init_cond.
void set_init_cond_manually(Vector& u, const InitialCondition& init_cond, const Grid& grid);

// Sets u (of length 2*grid.num_points()) from a frame Result-*.npy of the same grid (see read_frame).
// Without scaling the frame is copied as it is. With scaling every second point of the frame
// (a grid with half the resolution) is placed into the centre of the domain, the rest is the beta phase.
// Prints the reason to std::cerr and returns false on failure.
bool set_init_cond_from_file(Vector& u, const std::filesystem::path& filename, const Grid& grid, bool scaling);

// Makes the boundary points of the initial condition u satisfy the boundary conditions in params:
// Dirichlet sets the boundary value (phase beta = 0, concentration c_init_beta), Neumann copies the
// value of the neighbouring point. The time stepping does not write u on the boundary (Model gives
// Dirichlet points fu = 0 and Neumann points the fu of their neighbour), so the boundary conditions
// hold during the run only if they hold at the start. Call it after setting the initial condition.
void set_init_cond_boundary(Vector& u, const Grid& grid, const ModelParams& params);