#pragma once
#include <string>
#include <stdexcept>
#include "core/types.hpp"

enum class ICType {
    HyperbolicTangent,
    Circle, // sloped
    ConstantCircle,
    ConstantHalves,
    Stripe,
    TwoBumps,
    Star,
    ThreeBumps,
    PerpendicularStripes,
    Box,
    RandomBumps,
    WulffShape, // sloped
    Count 
};

enum class BCType {
    Dirichlet,
    Neumann
};

// The prescribed force term of AllenCahn (parameters.force_term).
enum class FTType {
    Constant,
    Radial
};

std::string BCType_to_string(BCType bc_type);
BCType string_to_BCType(const std::string& str);

std::string ICType_to_string(ICType ic_type);
ICType string_to_ICType(const std::string& str);

std::string FTType_to_string(FTType ft_type);
FTType string_to_FTType(const std::string& str);