#pragma once
#include "core/types.hpp"

// Constants shared by the whole code. The thermodynamic data of the alloy systems come from databases
// (databases/, read by io/TdbDatabase).
namespace constants {
    // M_PI is not standard C++ (glibc provides it, other headers may not); same literal as glibc's M_PI.
    inline constexpr Real pi = 3.14159265358979323846;

    // The gas constant [J/(mol K)] and calories to joules. They are double because only the host uses
    // them, to precompute the values that depend on the temperature (the mobilities of ZrNbChoud2013;
    // the Gibbs energies from a database use the database's R, TdbDatabase::R).
    // Do not use them in GPU code, it would compute in double: use the precomputed values there.
    inline constexpr double R = 8.31446261815324;
    inline constexpr double Cal_to_joule = 4.184;

    // The values of the phase field in the two phases.
    enum class Phase {
        alpha = 1,
        beta = 0
    };
}
