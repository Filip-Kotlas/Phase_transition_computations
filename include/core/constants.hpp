#pragma once
#include "core/types.hpp"

// Constants shared by the whole code. The data of the alloy systems are in models/thermodynamics/.
namespace constants {
    // M_PI is not standard C++ (glibc provides it, other headers may not); same literal as glibc's M_PI.
    inline constexpr Real pi = 3.14159265358979323846;

    // The gas constant [J/(mol K)] and calories to joules. They are double because only the host uses
    // them, to precompute the values that depend on the temperature (the thermodynamics structs).
    // Do not use them in GPU code, it would compute in double: use the precomputed values there.
    inline constexpr double R = 8.31446261815324;
    inline constexpr double Cal_to_joule = 4.184;

    // The values of the phase field in the two phases.
    enum class Phase {
        alpha = 1,
        beta = 0
    };
}
