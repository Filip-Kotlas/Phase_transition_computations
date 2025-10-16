#pragma once

#include <string>
#include <json.hpp>
#include <filesystem>
#include "constants.hpp"

class Parameters {
public:
    // data
    double initial_time{};
    double final_time{};
    Domain domain;
    int sizeX{};
    int sizeY{};
    int frame_num{};
    ICType init_condition{};
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

    // uložení do textového souboru ve formátu vhodném pro zkopírování do LateXu
    void save_for_latex(const std::filesystem::path& filename) const;

    // uloží přesnou kopii původního JSON
    void save_copy_of_config(const std::filesystem::path& original_path,
                             const std::filesystem::path& copy_path) const;

};