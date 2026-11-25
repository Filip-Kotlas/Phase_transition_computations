#include <cmath>
#include <string>

#include "models/ZrNbChoud2013.hpp"
#include "config/ConfigSection.hpp"
#include "core/constants.hpp"
#include "io/TdbDatabase.hpp"

namespace {
// The mobilities of Nb, functions of the temperature T [K], evaluated once, on the host, in double
// (the thermodynamic databases hold no diffusion data).
namespace from_T {
    using constants::R;
    using constants::Cal_to_joule;

    // the diffusivities of Nb in the two phases
    double D_Nb_alpha(double T) {
        return 6.6e-10*std::exp(-31500*Cal_to_joule/(R*T));
    }
    double D_Nb_beta(double T) {
        return 9e-10*std::pow(T/1136.0, 18.1)*std::exp(-(25100 + 35.5*(T - 1136))*Cal_to_joule/(R*T));
    }

    double M_Nb_alpha(double T) {
        return D_Nb_alpha(T) / R / T;
    }
    double M_Nb_beta(double T) {
        return D_Nb_beta(T) / R / T;
    }

    // Not used by the solver; kept for reference.
    // constexpr double delta = 5e-9;
    // constexpr double delta_fixed = 5e-9;
    // constexpr double molar_volume = 1.4060e-5;
    // constexpr double b = 3.23e-10;
    // constexpr double sigma = 0.3;
    // constexpr double l = 50*delta;
    // constexpr double w = 3*sigma/delta;
    // constexpr double c_min = 0.0001;
    // constexpr double c_max = 1 - c_min;
    // double D_eff(double T) {
    //     return std::pow(std::sqrt(D_Nb_alpha(T)) + std::sqrt(D_Nb_beta(T)), 2);
    // }
    // double M_phi(double T) {
    //     return 0.235*D_eff(T)*delta_fixed/delta*molar_volume/(b*b*R*T);
    // }
    // double M_phi_tilde(double T) {
    //     return M_phi(T)*l*l/(M_Nb_alpha(T)*molar_volume);
    // }
    // double epsilon() {
    //     return std::sqrt(6*sigma*delta);
    // }
    // double epsilon_tilde(double T) {
    //     return epsilon()/(l*std::sqrt(R*T/molar_volume));
    // }
    // double w_tilde(double T) {
    //     return w/(R*T);
    // }
}
}

ZrNbChoud2013::ZrNbChoud2013(const Parameters& parameters, const Params& own)
    : Model(parameters),
      own(own),
      M_Nb_alpha(Real(from_T::M_Nb_alpha(own.T))),
      M_Nb_beta(Real(from_T::M_Nb_beta(own.T)))
{
}

ZrNbChoud2013::Params ZrNbChoud2013::Params::load(const nlohmann::json& section)
{
    ConfigSection config(section, "parameters");
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
    params.d = config.get<Real>("d");
    params.T = config.get<Real>("T");
    const std::string database = config.get<std::string>("database");
    const std::string alpha_phase = config.get<std::string>("alpha_phase");
    const std::string beta_phase = config.get<std::string>("beta_phase");
    config.check_all_read();

    const TdbDatabase tdb = TdbDatabase::load(database);
    params.alpha_phase = tdb.binary_phase< SoluteSite::substitutional, 2 >(alpha_phase, "NB", "ZR", params.T);
    params.beta_phase = tdb.binary_phase< SoluteSite::substitutional, 2 >(beta_phase, "NB", "ZR", params.T);
    return params;
}

InitialConcentrations ZrNbChoud2013::Params::initial_concentrations() const
{
    return {c_init_alpha, c_init_beta};
}

std::vector<ParameterEntry> ZrNbChoud2013::Params::entries() const
{
    std::vector<ParameterEntry> result;
    result.push_back({"Alpha", "\\(\\alpha\\)", ParameterEntry::to_text(alpha)});
    result.push_back({"a", "\\(a\\)", ParameterEntry::to_text(par_a)});
    result.push_back({"b", "\\(b\\)", ParameterEntry::to_text(par_b)});
    result.push_back({"Ksi", "\\(\\xi\\)", ParameterEntry::to_text(ksi)});
    result.push_back({"A", "\\(A\\)", ParameterEntry::to_text(anisotropy.A)});
    result.push_back({"m", "\\(m\\)", ParameterEntry::to_text(anisotropy.m)});
    result.push_back({"theta_0", "\\(\\theta_0\\)", ParameterEntry::to_text(anisotropy.theta_0)});
    result.push_back({"d", "\\(d\\)", ParameterEntry::to_text(d)});
    result.push_back({"T", "\\(T\\)", ParameterEntry::to_text(T)});
    return result;
}

__cuda_callable__
void ZrNbChoud2013::set_rhs_at(const VectorView& u, VectorView& fu, Index i, Index j) const
{
    const Index offset = grid.num_points();
    fu[grid.index(i, j)] = get_rhs_phase_at(u, i, j);
    fu[offset + grid.index(i, j)] = get_rhs_concentration_at(u, i, j);
}

__cuda_callable__
Real ZrNbChoud2013::get_rhs_concentration_at(const VectorView& u, Index i, Index j) const
{
    return div_D_grad_concentration(u, i, j) + div_D_grad_phase(u, i, j);
}

__cuda_callable__
Real ZrNbChoud2013::F(const VectorView& u, Index i, Index j) const
{
    return own.alpha_phase.G(conc_at(u, i, j)) - own.beta_phase.G(conc_at(u, i, j));
}

__cuda_callable__
Real ZrNbChoud2013::div_D_grad_concentration(const VectorView& u, Index i, Index j) const
{
	Real coeff_plus_half = (get_conc_diff_coef(u, i + 1, j) + get_conc_diff_coef(u, i, j)) / 2;
	Real coeff_minus_half = (get_conc_diff_coef(u, i, j) + get_conc_diff_coef(u, i - 1, j)) / 2;
	Real x_direction = coeff_plus_half * (conc_at(u, i+1, j) - conc_at(u, i, j)) / grid.hx
				         - coeff_minus_half * (conc_at(u, i, j) - conc_at(u, i-1, j)) / grid.hx;

	coeff_plus_half = (get_conc_diff_coef(u, i, j+1) + get_conc_diff_coef(u, i, j)) / 2;
	coeff_minus_half = (get_conc_diff_coef(u, i, j) + get_conc_diff_coef(u, i, j-1)) / 2;
	Real y_direction = coeff_plus_half * (conc_at(u, i, j+1) - conc_at(u, i, j)) / grid.hy
                      	 - coeff_minus_half * (conc_at(u, i, j) - conc_at(u, i, j-1)) / grid.hy;

	return x_direction / grid.hx + y_direction / grid.hy;
}

__cuda_callable__
Real ZrNbChoud2013::div_D_grad_phase(const VectorView& u, Index i, Index j) const
{
	Real coeff_plus_half = (get_phas_diff_coef(u, i + 1, j) + get_phas_diff_coef(u, i, j)) / 2;
	Real coeff_minus_half = (get_phas_diff_coef(u, i, j) + get_phas_diff_coef(u, i - 1, j)) / 2;
	Real x_direction = coeff_plus_half * (phase_at(u, i+1, j) - phase_at(u, i, j)) / grid.hx
				         - coeff_minus_half * (phase_at(u, i, j) - phase_at(u, i-1, j)) / grid.hx;

	coeff_plus_half = (get_phas_diff_coef(u, i, j+1) + get_phas_diff_coef(u, i, j)) / 2;
	coeff_minus_half = (get_phas_diff_coef(u, i, j) + get_phas_diff_coef(u, i, j-1)) / 2;
	Real y_direction = coeff_plus_half * (phase_at(u, i, j+1) - phase_at(u, i, j)) / grid.hy
                      	 - coeff_minus_half * (phase_at(u, i, j) - phase_at(u, i, j-1)) / grid.hy;

	return x_direction / grid.hx + y_direction / grid.hy;
}

__cuda_callable__
Real ZrNbChoud2013::get_conc_diff_coef(const VectorView& u, Index i, Index j) const
{
    return own.par_b*own.d/Real(30.0)
            * conc_at(u, i, j)
            * (1 - conc_at(u, i, j))
            * M_Nb_beta
		    * pow(M_Nb_alpha/M_Nb_beta, polynom_p(u, i, j))
		    * sec_deriv_of_g_w_resp_to_c(u, i, j);
}

__cuda_callable__
Real ZrNbChoud2013::get_phas_diff_coef(const VectorView& u, Index i, Index j) const
{
   return own.par_b*own.d/Real(30.0)
          * conc_at(u, i, j)
		  * (1 - conc_at(u, i, j))
          * M_Nb_beta
		  * pow(M_Nb_alpha/M_Nb_beta, polynom_p(u, i, j))
		  * deriv_of_g_w_resp_to_c_and_p(u, i, j);
}

__cuda_callable__
Real ZrNbChoud2013::polynom_p(const VectorView& u, Index i, Index j) const
{
   return 6*pow(phase_at(u, i, j), Real(5)) - 15*pow(phase_at(u, i, j), Real(4)) + 10*pow(phase_at(u, i, j), Real(3));
}

__cuda_callable__
Real ZrNbChoud2013::der_polynom_p(const VectorView& u, Index i, Index j) const
{
   return 30*(pow(phase_at(u, i, j), Real(4)) - 2*pow(phase_at(u, i, j), Real(3)) + pow(phase_at(u, i, j), Real(2)));
}

__cuda_callable__
Real ZrNbChoud2013::sec_deriv_of_g_w_resp_to_c(const VectorView& u, Index i, Index j) const
{
   Real c = conc_at(u, i, j);
   return polynom_p(u, i, j)*own.alpha_phase.d2G_du2(c) + (1 - polynom_p(u, i, j))*own.beta_phase.d2G_du2(c);
}

__cuda_callable__
Real ZrNbChoud2013::deriv_of_g_w_resp_to_c_and_p(const VectorView& u, Index i, Index j) const
{
   Real c = conc_at(u, i, j);
   return der_polynom_p(u, i, j)*(own.alpha_phase.dG_du(c) - own.beta_phase.dG_du(c));
}

__cuda_callable__
Real ZrNbChoud2013::get_rhs_phase_at(const VectorView& u, Index i, Index j) const
{
    return Real(1.0)/own.alpha * (div_T0(u, i, j) + f_0(u, i , j) / own.ksi / own.ksi - own.par_b/own.ksi*grade_4_polynom(u, i, j)*F(u, i, j));
}

__cuda_callable__
Real ZrNbChoud2013::grad_p_1_forward(const VectorView& u, Index i, Index j) const
{
    return (phase_at(u, i+1, j) - phase_at(u, i, j))/grid.hx;
}

__cuda_callable__
Real ZrNbChoud2013::grad_p_2_forward(const VectorView& u, Index i, Index j) const
{
   return (phase_at(u, i, j+1) - phase_at(u, i, j))/grid.hy;
}

__cuda_callable__
Real ZrNbChoud2013::div_T0(const VectorView& u, Index i, Index j) const
{
    return (T0_1(grad_p_1_forward(u, i, j), grad_p_2_forward(u, i, j)) - T0_1(grad_p_1_forward(u, i-1, j), grad_p_2_forward(u, i-1, j)))/grid.hx
           + (T0_2(grad_p_1_forward(u, i, j), grad_p_2_forward(u, i, j)) - T0_2(grad_p_1_forward(u, i, j-1), grad_p_2_forward(u, i, j-1)))/grid.hy;
}

__cuda_callable__
Real ZrNbChoud2013::T0_1(const Real grad_p_1, const Real grad_p_2) const
{
    Real theta = atan2(-grad_p_2, -grad_p_1);
    const Anisotropy& a = own.anisotropy;
    return a.psi(theta) * a.psi(theta) * grad_p_1 - a.psi(theta) * a.der_psi(theta) * grad_p_2;
}

__cuda_callable__
Real ZrNbChoud2013::T0_2(const Real grad_p_1, const Real grad_p_2) const
{
    Real theta = atan2(-grad_p_2, -grad_p_1);
    const Anisotropy& a = own.anisotropy;
    return a.psi(theta) * a.psi(theta) * grad_p_2 + a.psi(theta) * a.der_psi(theta) * grad_p_1;
}

__cuda_callable__
Real ZrNbChoud2013::f_0(const VectorView& u, Index i, Index j) const
{
   return own.par_a*phase_at(u, i, j)*(1 - phase_at(u, i, j))*(phase_at(u, i, j) - Real(1.0)/Real(2.0));
}

__cuda_callable__
Real ZrNbChoud2013::grade_4_polynom(const VectorView& u, Index i, Index j) const
{
   return pow(phase_at(u, i, j), Real(2)) * pow(phase_at(u, i, j) - Real(1.0), Real(2));
}
