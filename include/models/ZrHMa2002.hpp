#pragma once

// The hydride precipitation in the Zr-H alloy after Ma et al. (2002) (model.type "Zr-H_Ma2002").
// So far only prepared: the parameters are read, written to parameters.txt/.tex and stored in the
// model, but the equations are not implemented yet (set_rhs_at gives zero for both fields, so the
// initial condition does not change). See Model.hpp for what a model provides.
//
// The model has no interface width ksi and no anisotropy, so it cannot use the initial conditions
// that need them (HyperbolicTangent, WulffShape): start in main.cu refuses them.
//
// Its config section (all keys required, all floating-point numbers):
//     "Zr-H_Ma2002": { "C1": ..., "C2": ..., "A1": ..., "A2": ..., "A3": ..., "A4": ...,
//                      "L_p": ..., "M": ... }

#include <vector>
#include "external/json.hpp"

#include "models/Model.hpp"

class ZrHMa2002 : public Model
{
   public:

    // The name of the model in the config: model.type and the section model.Zr-H_Ma2002.
    static constexpr const char* name = "Zr-H_Ma2002";

    // The physical parameters of this model (the config section model.Zr-H_Ma2002). Plain data,
    // copied to the GPU with the model.
    struct Params {
        Real C1;
        Real C2;
        Real A1;
        Real A2;
        Real A3;
        Real A4;
        Real L_p;
        Real M;

        // Reads the section strictly (ConfigSection.hpp); throws std::runtime_error on a wrong config.
        static Params load(const nlohmann::json& section);
        // The parameters for parameters.txt and parameters.tex.
        std::vector<ParameterEntry> entries() const;
    };

    ZrHMa2002(const Parameters& parameters, const Params& own) : Model(parameters), own(own) {}

    __cuda_callable__
    void set_rhs_at(const VectorView& u, VectorView& fu, Index i, Index j) const;

   private:

    const Params own;
};
