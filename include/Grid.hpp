#pragma once

// The uniform grid of sizeX x sizeY points on the rectangle domain, including the boundary points.
// The point (i, j), 0 <= i < sizeX, 0 <= j < sizeY, lies at (x(i), y(j)) and has the index
// index(i, j) = j*sizeX + i in each block of the state vector u (phase in [0, N), concentration in
// [N, 2N), N = points()). Plain data, built once from the config (Parameters::grid()) and copied into
// the GPU lambdas by value. hx and hy are derived from the other fields, so do not change those after
// the construction.

#include "types.hpp"

struct Domain {
    Real x_left;
    Real x_right;
    Real y_left;
    Real y_right;
};

struct Grid {
    Index sizeX = 0;
    Index sizeY = 0;
    Domain domain = {};
    Real hx = 0;  // (x_right - x_left) / (sizeX - 1)
    Real hy = 0;  // (y_right - y_left) / (sizeY - 1)

    Grid() = default;

    Grid(Index sizeX, Index sizeY, const Domain& domain)
    :  sizeX(sizeX),
       sizeY(sizeY),
       domain(domain),
       hx((domain.x_right - domain.x_left) / (sizeX - 1)),
       hy((domain.y_right - domain.y_left) / (sizeY - 1))
    {}

    // Number of grid points N (the length of one block of u).
    __cuda_callable__
    Index points() const { return sizeX * sizeY; }

    // Position of the point (i, j) in a block of u.
    __cuda_callable__
    Index index(Index i, Index j) const { return j * sizeX + i; }

    // Coordinates of the point (i, j).
    __cuda_callable__
    Real x(Index i) const { return domain.x_left + i * hx; }
    __cuda_callable__
    Real y(Index j) const { return domain.y_left + j * hy; }
};
