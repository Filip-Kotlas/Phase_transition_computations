#include <stdexcept>

#include "models/ZrHMa2002.hpp"
#include "config/ConfigSection.hpp"

ZrHMa2002::Params ZrHMa2002::Params::load(const nlohmann::json& section)
{
    ConfigSection config(section, "parameters");
    Params params;
    params.C1    = config.get<Real>("C1");
    params.C2    = config.get<Real>("C2");
    params.A1    = config.get<Real>("A1");
    params.A2    = config.get<Real>("A2");
    params.A3    = config.get<Real>("A3");
    params.alpha = config.get<Real>("alpha");
    params.beta  = config.get<Real>("beta");
    params.L_p   = config.get<Real>("L_p");
    params.M     = config.get<Real>("M");
    params.c_init_hydride = config.get<Real>("c_init_hydride");
    params.c_init_matrix  = config.get<Real>("c_init_matrix");
    config.check_all_read();
    // A negative coefficient turns a diffusion term into a backward one (an ill-posed problem).
    if (params.alpha < 0 || params.beta < 0 || params.L_p < 0 || params.M < 0)
        throw std::runtime_error(std::string("The parameters alpha, beta, L_p and M of model.") + name
                                 + " must not be negative.");
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
    result.push_back({"Alpha", "\\(\\alpha\\)", ParameterEntry::to_text(alpha)});
    result.push_back({"Beta", "\\(\\beta\\)", ParameterEntry::to_text(beta)});
    result.push_back({"L_p", "\\(L_p\\)", ParameterEntry::to_text(L_p)});
    result.push_back({"M", "\\(M\\)", ParameterEntry::to_text(M)});
    result.push_back({"Initial concentration of the hydride", "\\(C\\) init. hydride", ParameterEntry::to_text(c_init_hydride)});
    result.push_back({"Initial concentration of the matrix", "\\(C\\) init. matrix", ParameterEntry::to_text(c_init_matrix)});
    return result;
}

InitialConcentrations ZrHMa2002::Params::initial_concentrations() const
{
    return {c_init_hydride, c_init_matrix};
}

void ZrHMa2002::check_parameters(const Parameters& parameters)
{
    // The boundary of the Cahn-Hilliard equation is the no-flux pair dC/dn = 0, dmu/dn = 0 (the clamp in
    // chemical_potential_at); a fixed concentration on the boundary has not been thought through.
    if (parameters.bound_cond.conc_x != BCType::Neumann || parameters.bound_cond.conc_y != BCType::Neumann)
        throw std::runtime_error(std::string("The model ") + name
                                 + " needs Neumann boundary conditions for the concentration (bound_cond.x_conc and y_conc).");
}

__cuda_callable__
void ZrHMa2002::set_rhs_at(const VectorView& u, VectorView& fu, Index i, Index j) const
{
    const Index offset = grid.num_points();
    fu[grid.index(i, j)] = get_rhs_phase_at(u, i, j);
    fu[offset + grid.index(i, j)] = get_rhs_concentration_at(u, i, j);
}

// d phi/dt = L_p (alpha lap phi - df/dphi),  df/dphi = A2 (C - C2) phi - A3 phi^3 + A3 phi^5
__cuda_callable__
Real ZrHMa2002::get_rhs_phase_at(const VectorView& u, Index i, Index j) const
{
    const Real p = phase_at(u, i, j);
    const Real c = conc_at(u, i, j);
    const Real p2 = p * p;
    const Real df_dp = p * (own.A2 * (c - own.C2) - own.A3 * p2 * (Real(1) - p2));
    return own.L_p * (own.alpha * laplace_phase_at(u, i, j) - df_dp);
}

// dC/dt = M lap mu
__cuda_callable__
Real ZrHMa2002::get_rhs_concentration_at(const VectorView& u, Index i, Index j) const
{
    const Real mu = chemical_potential_at(u, i, j);
    const Real hx2 = grid.hx * grid.hx;
    const Real hy2 = grid.hy * grid.hy;
    return own.M * ((chemical_potential_at(u, i + 1, j) - Real(2) * mu + chemical_potential_at(u, i - 1, j)) / hx2
                    + (chemical_potential_at(u, i, j + 1) - Real(2) * mu + chemical_potential_at(u, i, j - 1)) / hy2);
}

// mu = A1 (C - C1) + A2/2 phi^2 - beta lap C
__cuda_callable__
Real ZrHMa2002::chemical_potential_at(const VectorView& u, Index i, Index j) const
{
    // The no-flux condition dmu/dn = 0: mu on the boundary equals mu at the neighbouring interior point.
    // (Only the edge points are ever asked for, never the corners or anything outside the grid.)
    if (i == 0) i = 1;
    else if (i == grid.sizeX - 1) i = grid.sizeX - 2;
    if (j == 0) j = 1;
    else if (j == grid.sizeY - 1) j = grid.sizeY - 2;

    const Real p = phase_at(u, i, j);
    const Real c = conc_at(u, i, j);
    return own.A1 * (c - own.C1) + own.A2 / Real(2) * p * p - own.beta * laplace_conc_at(u, i, j);
}

__cuda_callable__
Real ZrHMa2002::laplace_phase_at(const VectorView& u, Index i, Index j) const
{
    const Real p = phase_at(u, i, j);
    return (phase_at(u, i + 1, j) - Real(2) * p + phase_at(u, i - 1, j)) / (grid.hx * grid.hx)
           + (phase_at(u, i, j + 1) - Real(2) * p + phase_at(u, i, j - 1)) / (grid.hy * grid.hy);
}

__cuda_callable__
Real ZrHMa2002::laplace_conc_at(const VectorView& u, Index i, Index j) const
{
    const Real c = conc_at(u, i, j);
    return (conc_at(u, i + 1, j) - Real(2) * c + conc_at(u, i - 1, j)) / (grid.hx * grid.hx)
           + (conc_at(u, i, j + 1) - Real(2) * c + conc_at(u, i, j - 1)) / (grid.hy * grid.hy);
}
