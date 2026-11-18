#include <cmath>
#include <algorithm>

#include "AllenCahnModel.hpp"
#include "ConfigSection.hpp"

AllenCahnModel::Params AllenCahnModel::Params::load(const nlohmann::json& section)
{
    ConfigSection config(section, std::string("model.") + name);
    Params params;
    // the phase equation
    params.alpha = config.get<Real>("alpha");
    params.par_a = config.get<Real>("a");
    params.par_b = config.get<Real>("b");
    params.ksi   = config.get<Real>("ksi");
    // m has to be an integer, otherwise psi is not 2*pi-periodic (get<int> refuses 4.5).
    params.anisotropy.A       = config.get<Real>("A");
    params.anisotropy.m       = config.get<int>("m");
    params.anisotropy.theta_0 = config.get<Real>("theta_0");
    const int m = params.anisotropy.m;
    if (m < 2)
        throw std::runtime_error("The anisotropy parameter m has to be at least 2.");
    if (params.anisotropy.A > 1.0/(m*m - 1))
        throw std::runtime_error("Parameters A and m do not fulfill the convexity condition.");
    const std::string force_term = config.get<std::string>("force_term");
    if (force_term != "constant" && force_term != "radial")
        throw std::runtime_error("The config value model.allen_cahn.force_term has to be \"constant\" or \"radial\", it is \""
                                 + force_term + "\".");
    params.force_term = string_to_FTType(force_term);
    // The size only means something for the constant force term; check_all_read rejects it otherwise.
    params.force_size = (params.force_term == FTType::Constant) ? config.get<Real>("force_size") : Real(0);
    config.check_all_read();
    return params;
}

std::vector<ParameterEntry> AllenCahnModel::Params::entries() const
{
    std::vector<ParameterEntry> result;
    result.push_back({"Alpha", "\\(\\alpha\\)", ParameterEntry::to_text(alpha)});
    result.push_back({"a", "\\(a\\)", ParameterEntry::to_text(par_a)});
    result.push_back({"b", "\\(b\\)", ParameterEntry::to_text(par_b)});
    result.push_back({"Ksi", "\\(\\xi\\)", ParameterEntry::to_text(ksi)});
    result.push_back({"A", "\\(A\\)", ParameterEntry::to_text(anisotropy.A)});
    result.push_back({"m", "\\(m\\)", ParameterEntry::to_text(anisotropy.m)});
    result.push_back({"theta_0", "\\(\\theta_0\\)", ParameterEntry::to_text(anisotropy.theta_0)});
    result.push_back({"Force term", "force term", FTType_to_string(force_term)});
    if (force_term == FTType::Constant)
        result.push_back({"Force term size", "force term size", ParameterEntry::to_text(force_size)});
    return result;
}

__cuda_callable__
void AllenCahnModel::set_rhs_at(const VectorView& u, VectorView& fu, Index i, Index j) const
{
    const Index offset = grid.num_points();
    fu[grid.index(i, j)] = get_rhs_phase_at(u, i, j, F(i, j));
    fu[offset + grid.index(i, j)] = 0;
}

__cuda_callable__
Real AllenCahnModel::F(Index i, Index j) const
{
    if (own.force_term == FTType::Constant)
        return own.force_size;
    else if (own.force_term == FTType::Radial)
    {
        Real mid_x = (grid.domain.x_right - grid.domain.x_left)/2;
        Real mid_y = (grid.domain.y_right - grid.domain.y_left)/2;
        Real r = sqrt(pow(i*grid.hx - mid_x, Real(2)) + pow(j*grid.hy - mid_y, Real(2)));
        return -2/std::max(r, Real(0.1));
    }
    else
        return 0;
}

__cuda_callable__
Real AllenCahnModel::get_rhs_phase_at(const VectorView& u, Index i, Index j, Real F) const
{
    return Real(1.0)/own.alpha * (div_T0(u, i, j) + f_0(u, i , j) / own.ksi / own.ksi - own.par_b/own.ksi*grade_4_polynom(u, i, j)*F);
}

__cuda_callable__
Real AllenCahnModel::grad_p_1_forward(const VectorView& u, Index i, Index j) const
{
    return (phase_at(u, i+1, j) - phase_at(u, i, j))/grid.hx;
}

__cuda_callable__
Real AllenCahnModel::grad_p_2_forward(const VectorView& u, Index i, Index j) const
{
   return (phase_at(u, i, j+1) - phase_at(u, i, j))/grid.hy;
}

__cuda_callable__
Real AllenCahnModel::div_T0(const VectorView& u, Index i, Index j) const
{
    return (T0_1(grad_p_1_forward(u, i, j), grad_p_2_forward(u, i, j)) - T0_1(grad_p_1_forward(u, i-1, j), grad_p_2_forward(u, i-1, j)))/grid.hx
           + (T0_2(grad_p_1_forward(u, i, j), grad_p_2_forward(u, i, j)) - T0_2(grad_p_1_forward(u, i, j-1), grad_p_2_forward(u, i, j-1)))/grid.hy;
}

__cuda_callable__
Real AllenCahnModel::T0_1(const Real grad_p_1, const Real grad_p_2) const
{
    Real theta = atan2(-grad_p_2, -grad_p_1);
    const Anisotropy& a = own.anisotropy;
    return a.psi(theta) * a.psi(theta) * grad_p_1 - a.psi(theta) * a.der_psi(theta) * grad_p_2;
}

__cuda_callable__
Real AllenCahnModel::T0_2(const Real grad_p_1, const Real grad_p_2) const
{
    Real theta = atan2(-grad_p_2, -grad_p_1);
    const Anisotropy& a = own.anisotropy;
    return a.psi(theta) * a.psi(theta) * grad_p_2 + a.psi(theta) * a.der_psi(theta) * grad_p_1;
}

__cuda_callable__
Real AllenCahnModel::f_0(const VectorView& u, Index i, Index j) const
{
   return own.par_a*phase_at(u, i, j)*(1 - phase_at(u, i, j))*(phase_at(u, i, j) - Real(1.0)/Real(2.0));
}

__cuda_callable__
Real AllenCahnModel::grade_4_polynom(const VectorView& u, Index i, Index j) const
{
   return pow(phase_at(u, i, j), Real(2)) * pow(phase_at(u, i, j) - Real(1.0), Real(2));
}
