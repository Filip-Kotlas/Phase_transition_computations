#include "enums.hpp"

std::string BCType_to_string(BCType bc_type) {
    switch (bc_type) {
        case BCType::Dirichlet:
            return "dirichlet";
        case BCType::Neumann:
            return "neumann";
        default:
            return "unknown";
    }
}

BCType string_to_BCType(const std::string& str) {
    if (str == "dirichlet") {
        return BCType::Dirichlet;
    } else if (str == "neumann") {
        return BCType::Neumann;
    } else {
        throw std::runtime_error("Unknown boundary condition type: " + str);
    }
}

std::string ICType_to_string(ICType ic_type) {
    switch (ic_type) {
        case ICType::HyperbolicTangent:
            return "HyperbolicTangent";
        case ICType::LinearByParts:
            return "LinearByParts";
        case ICType::ConstantCircle:
            return "ConstantCircle";
        case ICType::ConstantHalves:
            return "ConstantHalves";
        case ICType::Stripe:
            return "Stripe";
        case ICType::TwoBumps:
            return "TwoBumps";
        case ICType::Star:
            return "Star";
        case ICType::ThreeBumps:
            return "ThreeBumps";
        case ICType::PerpendicularStripes:
            return "PerpendicularStripes";
        case ICType::Box:
            return "Box";
        case ICType::RandomBumps:
            return "RandomBumps";
        case ICType::WulffShape:
            return "WulffShape";
        default:
            return "unknown";
    }
}

ICType string_to_ICType(const std::string& str) {
    if (str == "HyperbolicTangent") {
        return ICType::HyperbolicTangent;
    } else if (str == "LinearByParts") {
        return ICType::LinearByParts;
    } else if (str == "ConstantCircle") {
        return ICType::ConstantCircle;
    } else if (str == "ConstantHalves") {
        return ICType::ConstantHalves;
    } else if (str == "Stripe") {
        return ICType::Stripe;
    } else if (str == "TwoBumps") {
        return ICType::TwoBumps;
    } else if (str == "Star") {
        return ICType::Star;
    } else if (str == "ThreeBumps") {
        return ICType::ThreeBumps;
    } else if (str == "PerpendicularStripes") {
        return ICType::PerpendicularStripes;
    } else if (str == "Box") {
        return ICType::Box;
    } else if (str == "RandomBumps") {
        return ICType::RandomBumps;
    } else if (str == "WulffShape") {
        return ICType::WulffShape;
    } else {
        throw std::runtime_error("Unknown initial condition type: " + str);
    }
}

std::string FTType_to_string(FTType ft_type) {
    switch (ft_type) {
        case FTType::Constant:
            return "constant";
        case FTType::Zirconium:
            return "zirconium";
        case FTType::Radial:
            return "radial";
        default:
            return "unknown";
    }
}

FTType string_to_FTType(const std::string& str) {
    if (str == "constant") {
        return FTType::Constant;
    } else if (str == "zirconium") {
        return FTType::Zirconium;
    } else if (str == "radial") {
        return FTType::Radial;
    } else {
        throw std::runtime_error("Unknown force term type: " + str);
    }
}