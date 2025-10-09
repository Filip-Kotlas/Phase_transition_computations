#pragma once

#include <string>
#include <json.hpp>
#include <filesystem>

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

enum class InitialCondition {
    HyperbolicTangent,
    LinearByParts,
    ConstantCircle,
    ConstantHalves,
    Stripe,
    TwoBumps,
    Star,
    FourierX,
    FourierY
};

class Parameters {
public:
    // data
    std::string type;
    double initial_time{};
    double final_time{};
    Domain domain;
    int sizeX{};
    int sizeY{};
    int frame_num{};
    InitialCondition init_condition{InitialCondition::HyperbolicTangent};
    bool init_cond_from_file{};
    std::string init_cond_file_path;
    double timeStep{};
    double integrationTimeStep{};
    double alpha{};
    double beta{};
    double par_a{};
    double par_b{};
    double par_d{};
    double T{};
    double ksi{};
    MODEL model{MODEL::MODEL_1};

    // --- metody ---
    Parameters() = default;

    // načtení z JSON souboru
    static Parameters load(const std::filesystem::path& filename);

    // uložení do textového souboru (pro info/parameters.txt)
    void save_human_readable(const std::filesystem::path& filename) const;

    // uloží přesnou kopii původního JSON
    void save_copy_of_config(const std::filesystem::path& original_path,
                             const std::filesystem::path& copy_path) const;

};
