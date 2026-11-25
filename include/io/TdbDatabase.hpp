#pragma once

// A thermodynamic database in the TDB format (Thermo-Calc's database format, used by CALPHAD
// assessments), read on the host. It reads a strict subset: the statements ELEMENT, PHASE, CONSTITUENT,
// FUNCTION, PARAMETER (G and L), TYPE_DEFINITION and TEMPERATURE_LIMITS; the bookkeeping statements
// (SPECIES, DEFINE_SYSTEM_DEFAULT, DEFAULT_COMMAND, ASSESSED_SYSTEMS, LIST_OF_REFERENCES, ...) are
// skipped. Keywords may be abbreviated as in Thermo-Calc (PAR, CONST, TYPE-DEF, TEMP-LIM, ...); text after
// $ is a comment. Anything else is an error, and so is anything in a requested phase that is not
// understood (other parameter types such as the magnetic TC and BMAGN, type definitions other than SEQ and
// MAGNETIC, wildcards), so a database is never read differently from how Thermo-Calc would read it.
//
// The expressions of FUNCTION and PARAMETER are piecewise in T (the ranges low <= T < high, the last one
// including its upper limit; a ,, limit is the database's TEMPERATURE_LIMITS) and are built from numbers,
// T, P, + - * / **, LN, LOG (natural, as in Thermo-Calc), EXP, parentheses and the names of FUNCTIONs.
// P is 1e5 Pa (it appears only in the gas phases).
//
// binary_phase evaluates everything at a temperature, in double, and returns the plain data the GPU needs
// (BinaryPhase.hpp: the variable u, the reference states and the formula are described there). Checked
// against pycalphad 0.11 on the Zr-Nb and Zr-H databases of B. Hallstedt's collection (databases/).

#include <filesystem>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include "models/BinaryPhase.hpp"

class TdbDatabase
{
public:
    // The gas constant of the Gibbs energies of mixing, Thermo-Calc's value (the assessments were made
    // with it; constants::R = 8.31446... is the current CODATA value, pycalphad uses 8.3145).
    static constexpr double R = 8.31451;
    // The pressure [Pa] in the expressions.
    static constexpr double P = 1e5;

    // Reads the database; throws std::runtime_error when it cannot be read or uses unsupported syntax.
    static TdbDatabase load(const std::filesystem::path& file);

    // The phase `phase` of the binary system of the elements `solute` and `host` at the temperature T [K]
    // (element and phase names as in the database, case-insensitive). The phase has to have exactly one
    // sublattice where two constituents mix: (solute, host) with only vacancies on the others for
    // SoluteSite::substitutional, (solute, VA) with the host alone on the only other sublattice for
    // SoluteSite::interstitial; its Redlich-Kister terms are at most L_(terms-1). Throws std::runtime_error
    // otherwise, or when an end member is missing or T is outside the range of an expression.
    template< SoluteSite site, int terms >
    BinaryPhase< site, terms > binary_phase(const std::string& phase, const std::string& solute,
                                            const std::string& host, double T) const
    {
        const PhaseValues values = phase_values(phase, solute, host, T, site, terms);
        BinaryPhase< site, terms > result;
        result.k = Real(values.k);
        result.G_solute = Real(values.G_solute);
        result.G_other = Real(values.G_other);
        result.mixing_RT = Real(values.mixing_RT);
        for (int n = 0; n < terms; ++n)
            result.L[n] = Real(values.L[n]);
        return result;
    }

    // The value of the FUNCTION `name` at T (for checks).
    double function_value(const std::string& name, double T) const;

    const std::filesystem::path& file() const { return path; }

    // The parsed expressions (public only for the implementation in the .cu).
    struct Node;
    struct Range {
        double low, high;
        std::shared_ptr<const Node> expression;
    };
    using Piecewise = std::vector<Range>;

private:
    struct Phase {
        std::string type_codes;                       // the letters after the name, e.g. "%A"
        std::vector<double> sites;                    // per sublattice
        std::vector<std::vector<std::string>> constituents;  // per sublattice (empty until CONSTITUENT)
    };
    struct Parameter {
        std::string type;                             // G, L, TC, BMAGN, ...
        std::string phase;
        std::vector<std::vector<std::string>> constituents;  // per sublattice
        int order;
        Piecewise value;
        std::string text;                             // the statement, for error messages
    };

    // What binary_phase returns, in double (L has `terms` entries, the missing ones 0). Checks everything,
    // including the kind `site` and the highest order terms - 1.
    struct PhaseValues {
        double k, G_solute, G_other, mixing_RT;
        std::vector<double> L;
    };
    PhaseValues phase_values(const std::string& phase, const std::string& solute, const std::string& host,
                             double T, SoluteSite site, int terms) const;

    void parse_statement(const std::string& statement);
    Piecewise parse_piecewise(const std::string& text, const std::string& what) const;
    double evaluate(const Piecewise& piecewise, double T, const std::string& what, int depth) const;
    double evaluate(const Node& node, double T, int depth) const;
    double reference_energy(const std::string& element, double T) const;
    [[noreturn]] void fail(const std::string& message) const;

    std::filesystem::path path;
    double T_low = 298.15, T_high = 6000;   // TEMPERATURE_LIMITS, the defaults of Thermo-Calc
    std::set<std::string> elements;
    std::map<std::string, Phase> phases;
    std::map<std::string, std::string> type_definitions;   // letter -> the statement
    std::map<std::string, Piecewise> functions;
    std::vector<Parameter> parameters;
};
