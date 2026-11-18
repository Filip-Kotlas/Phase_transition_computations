#pragma once
#include <string>
#include <stdexcept>
#include "types.hpp"

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
    Radial
};

// The model (the class in main.cu's run<ModelT>) of a run; given by the force term type,
// see Parameters::model_type().
enum class ModelType {
    AllenCahn,  // AllenCahnModel
    Zirconium   // ZirconiumModel
};

std::string BCType_to_string(BCType bc_type);
BCType string_to_BCType(const std::string& str);

std::string ICType_to_string(ICType ic_type);
ICType string_to_ICType(const std::string& str);

std::string FTType_to_string(FTType ft_type);
FTType string_to_FTType(const std::string& str);