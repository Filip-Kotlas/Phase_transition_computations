#pragma once
#include <string>
#include <stdexcept>

struct Domain {
    double x_left;
    double x_right;
    double y_left;
    double y_right;
};

enum class MODEL {
    MODEL_1 = 1,
    MODEL_2 = 2,
    MODEL_3 = 3,
    MODEL_4 = 4
};

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

enum class FTType {
    Constant,
    Zirconium,
    Radial,
    Reality,
};

std::string BCType_to_string(BCType bc_type);
BCType string_to_BCType(const std::string& str);

std::string ICType_to_string(ICType ic_type);
ICType string_to_ICType(const std::string& str);

std::string FTType_to_string(FTType ft_type);
FTType string_to_FTType(const std::string& str);