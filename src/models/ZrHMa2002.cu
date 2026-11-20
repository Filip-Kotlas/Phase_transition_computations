#include "models/ZrHMa2002.hpp"
#include "config/ConfigSection.hpp"

ZrHMa2002::Params ZrHMa2002::Params::load(const nlohmann::json& section)
{
    ConfigSection config(section, std::string("model.") + name);
    Params params;
    params.C1  = config.get<Real>("C1");
    params.C2  = config.get<Real>("C2");
    params.A1  = config.get<Real>("A1");
    params.A2  = config.get<Real>("A2");
    params.A3  = config.get<Real>("A3");
    params.A4  = config.get<Real>("A4");
    params.L_p = config.get<Real>("L_p");
    params.M   = config.get<Real>("M");
    config.check_all_read();
    return params;
}

std::vector<ParameterEntry> ZrHMa2002::Params::entries() const
{
    std::vector<ParameterEntry> result;
    result.push_back({"C1", "\\(C_1\\)", ParameterEntry::to_text(C1)});
    result.push_back({"C2", "\\(C_2\\)", ParameterEntry::to_text(C2)});
    result.push_back({"A1", "\\(A_1\\)", ParameterEntry::to_text(A1)});
    result.push_back({"A2", "\\(A_2\\)", ParameterEntry::to_text(A2)});
    result.push_back({"A3", "\\(A_3\\)", ParameterEntry::to_text(A3)});
    result.push_back({"A4", "\\(A_4\\)", ParameterEntry::to_text(A4)});
    result.push_back({"L_p", "\\(L_p\\)", ParameterEntry::to_text(L_p)});
    result.push_back({"M", "\\(M\\)", ParameterEntry::to_text(M)});
    return result;
}

// TODO: the equations of the model. Until then both fields stay as they are.
__cuda_callable__
void ZrHMa2002::set_rhs_at(const VectorView& u, VectorView& fu, Index i, Index j) const
{
    const Index offset = grid.num_points();
    fu[grid.index(i, j)] = 0;
    fu[offset + grid.index(i, j)] = 0;
}
