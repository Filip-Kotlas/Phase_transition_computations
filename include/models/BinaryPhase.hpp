#pragma once

// The molar Gibbs energy of one phase of a binary system as a function of the concentration u of the
// solute, built on the host from a CALPHAD database (TdbDatabase::binary_phase) at a fixed temperature
// and evaluated on the GPU. Plain data: a model holds it, so it goes to the GPU with the model.
//
// The variable u is the u-fraction of the solute (Hillert): the amount of solute per atom of the host
// lattice, the same in every phase.
//   - substitutional solute, a phase (S,H)_a(Va)_c (Zr-Nb: (Nb,Zr)1(Va)3): u = n_S / (n_S + n_H), the mole
//     fraction of the solute; the site fraction of the solute on the mixing sublattice is y = u;
//   - interstitial solute, a phase (H)_a(S,Va)_c (Zr-H: (Zr)1(H,Va)2): u = n_S / n_H (the H/Zr ratio);
//     y = k*u with k = a/c.
// G is per mole of host-lattice atoms (per mole of atoms for a substitutional solute) [J/mol], relative to
// the reference states of the elements (each element's GHSER function at the same T is subtracted:
// the same linear function of u for all phases, so no equilibrium, driving force G_1(u) - G_2(u) or
// gradient of dG/du changes, and the float values stay small).
//
// In the compound energy formalism with a Redlich-Kister excess term, per formula unit:
//   G = y G_S + (1-y) G_O + c_mix R T (y ln y + (1-y) ln(1-y)) + y (1-y) sum_k L_k (y_S - y_O)^k,
// with G_S and G_O the end members with the solute and with the other constituent of the mixing
// sublattice (the host element, or the vacancy), and c_mix the sites of that sublattice. The members below
// hold the end members, R T c_mix and the L_k already divided by the host atoms per formula unit.
// The derivatives: G'(u) = y' g'(y), G''(u) = y'^2 g''(y) (y' = 1 or k), where, with d = 2y - 1,
//   g'(y)  = G_S - G_O + m (ln y - ln(1-y)) - d S(d) + 2 y(1-y) S'(d)
//   g''(y) = m / (y(1-y)) - 2 S(d) - 4 d S'(d) + 4 y(1-y) S''(d)
// (m = mixing_RT, S(d) = sum_k L_k d^k). For L_0 and L_1 alone this is
//   g'  = G_S - G_O + m (ln y - ln(1-y)) + (1-2y) L_0 + (6y - 6y^2 - 1) L_1,
//   g'' = m / (y(1-y)) - 2 L_0 + (6 - 12y) L_1.
//
// Two template parameters, both chosen by the model, make the GPU code as cheap as the hand-written one
// (bench vs 033d1ae on zirconium_200/400, where the first generic version ran 0.78/0.70x):
//   - site: for a substitutional solute y = u is known at compile time, so the phases of a model take the
//     logarithms of the same number and the compiler computes them once for all phases (with y = k*u and
//     k known only at run time it computed them per phase);
//   - terms: the number of Redlich-Kister terms L_0 ... L_(terms-1), so the compiler knows the higher ones
//     do not exist (unused terms that are 0 at run time still cost an FMA, a multiply-add instruction, on
//     every call; S'' and its term in g'' vanish entirely for terms <= 2). The reader refuses a phase of
//     the database with a higher order, so no term is dropped silently; a phase with fewer gets zeros.
// The code is in this header because a template cannot be compiled once for all its parameters.

#include <cmath>
#include "core/types.hpp"

// Where the solute sits: on the sublattice of the host element, or on an interstitial sublattice with
// vacancies.
enum class SoluteSite { substitutional, interstitial };

template< SoluteSite site, int terms >
struct BinaryPhase
{
    static_assert(terms >= 1, "a phase needs at least L_0 (it may be 0)");

    Real k;                  // y = k * u for an interstitial solute (not used for a substitutional one)
    Real G_solute;           // the end member with the solute on the mixing sublattice
    Real G_other;            // the end member with the other constituent there
    Real mixing_RT;          // R T * (sites of the mixing sublattice) / (host atoms per formula unit)
    Real L[terms];           // the Redlich-Kister parameters of (solute, other); missing ones are 0

    // G(u) and its first and second derivative with respect to u, for 0 < y < 1
    // (the y ln y terms give NaN at y = 0 and 1 in floating point).
    __cuda_callable__ Real G(Real u) const
    {
        const Real y = site_fraction(u);
        return y*G_solute
               + (1-y)*G_other
               + mixing_RT*(y*log(y) + (1-y)*log(1-y))
               + y*(1-y)*excess_sum(2*y - 1);
    }

    __cuda_callable__ Real dG_du(Real u) const
    {
        const Real y = site_fraction(u);
        const Real d = 2*y - 1;
        return dy_du()*(G_solute - G_other
                        + mixing_RT*(log(y) - log(1-y))
                        - d*excess_sum(d)
                        + 2*y*(1-y)*excess_sum_d1(d));
    }

    __cuda_callable__ Real d2G_du2(Real u) const
    {
        const Real y = site_fraction(u);
        const Real d = 2*y - 1;
        // mixing_RT times the reciprocal, not divided by y(1-y): the reciprocal is then the same expression
        // in all phases with the same y, so the compiler divides once (a float division is a subroutine
        // call on the GPU; dividing per phase cost a few per cent in bench TdbBench3).
        Real result = mixing_RT*(Real(1)/(y*(1-y)))
                      - 2*excess_sum(d)
                      - 4*d*excess_sum_d1(d);
        if constexpr (terms > 2)
            result += 4*y*(1-y)*excess_sum_d2(d);
        return dy_du()*dy_du()*result;
    }

private:
    // y(u) and dy/du
    __cuda_callable__ Real site_fraction(Real u) const
    {
        if constexpr (site == SoluteSite::substitutional)
            return u;
        else
            return k*u;
    }

    __cuda_callable__ Real dy_du() const
    {
        if constexpr (site == SoluteSite::substitutional)
            return 1;
        else
            return k;
    }

    // The Redlich-Kister sum S(d) = sum_k L_k d^k and its first and second derivative, by Horner's scheme.
    // The loops have a bound known at compile time, so they are unrolled and L stays in registers (a loop up
    // to a count known only at run time put L on the stack: bench TdbBench, 0.70-0.78x).
    __cuda_callable__ Real excess_sum(Real d) const
    {
        Real sum = 0;
        #pragma unroll
        for (int n = terms - 1; n >= 0; --n)
            sum = sum*d + L[n];
        return sum;
    }

    __cuda_callable__ Real excess_sum_d1(Real d) const
    {
        Real sum = 0;
        #pragma unroll
        for (int n = terms - 1; n >= 1; --n)
            sum = sum*d + n*L[n];
        return sum;
    }

    __cuda_callable__ Real excess_sum_d2(Real d) const
    {
        Real sum = 0;
        #pragma unroll
        for (int n = terms - 1; n >= 2; --n)
            sum = sum*d + n*(n-1)*L[n];
        return sum;
    }
};

template< int terms >
using SubstitutionalPhase = BinaryPhase< SoluteSite::substitutional, terms >;
template< int terms >
using InterstitialPhase = BinaryPhase< SoluteSite::interstitial, terms >;
