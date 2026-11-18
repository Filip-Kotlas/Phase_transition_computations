#include <cmath>

#include "ZirconiumModel.hpp"
#include "constants.hpp"
#include "ConfigSection.hpp"

ZirconiumModel::Params ZirconiumModel::Params::load(const nlohmann::json& section)
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
    params.d = config.get<Real>("d");
    params.T = config.get<Real>("T");
    config.check_all_read();
    return params;
}

std::vector<ParameterEntry> ZirconiumModel::Params::entries() const
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
void ZirconiumModel::set_rhs_at(const VectorView& u, VectorView& fu, Index i, Index j) const
{
    const Index offset = grid.num_points();
    fu[grid.index(i, j)] = get_rhs_phase_at(u, i, j, F(u, i, j));
    fu[offset + grid.index(i, j)] = get_rhs_concentration_at(u, i, j);
}

__cuda_callable__
Real ZirconiumModel::get_rhs_concentration_at(const VectorView& u, Index i, Index j) const
{
    return div_D_grad_concentration(u, i, j) + div_D_grad_phase(u, i, j);
}

__cuda_callable__
Real ZirconiumModel::F(const VectorView& u, Index i, Index j) const
{
    return constants::G_m_alpha(conc_at(u, i, j), own.T) - constants::G_m_beta(conc_at(u, i, j), own.T);
}

__cuda_callable__
Real ZirconiumModel::div_D_grad_concentration(const VectorView& u, Index i, Index j) const
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
Real ZirconiumModel::div_D_grad_phase(const VectorView& u, Index i, Index j) const
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
Real ZirconiumModel::get_conc_diff_coef(const VectorView& u, Index i, Index j) const
{
    return own.par_b*own.d/Real(30.0)
            * conc_at(u, i, j)
            * (1 - conc_at(u, i, j))
            * constants::M_Nb_beta(own.T)
		    * pow(constants::M_Nb_alpha(own.T)/constants::M_Nb_beta(own.T), polynom_p(u, i, j))
		    * sec_deriv_of_g_w_resp_to_c(u, i, j);
}

__cuda_callable__
Real ZirconiumModel::get_phas_diff_coef(const VectorView& u, Index i, Index j) const
{
   return own.par_b*own.d/Real(30.0)
          * conc_at(u, i, j)
		  * (1 - conc_at(u, i, j))
          * constants::M_Nb_beta(own.T)
		  * pow(constants::M_Nb_alpha(own.T)/constants::M_Nb_beta(own.T), polynom_p(u, i, j))
		  * deriv_of_g_w_resp_to_c_and_p(u, i, j);
}

__cuda_callable__
Real ZirconiumModel::polynom_p(const VectorView& u, Index i, Index j) const
{
   return 6*pow(phase_at(u, i, j), Real(5)) - 15*pow(phase_at(u, i, j), Real(4)) + 10*pow(phase_at(u, i, j), Real(3));
}

__cuda_callable__
Real ZirconiumModel::der_polynom_p(const VectorView& u, Index i, Index j) const
{
   return 30*(pow(phase_at(u, i, j), Real(4)) - 2*pow(phase_at(u, i, j), Real(3)) + pow(phase_at(u, i, j), Real(2)));
}

__cuda_callable__
Real ZirconiumModel::sec_deriv_of_g_w_resp_to_c(const VectorView& u, Index i, Index j) const
{
   Real c = conc_at(u, i, j);
   Real d2_G_alpha_wrt_c = constants::R * own.T / (c*(1-c))
   							 - 2 * constants::L_0_alpha;
   Real d2_G_beta_wrt_c = constants::R * own.T / (c*(1-c))
   							 - 2 * constants::L_0_beta(own.T)
							 + (6 - 12*c) * constants::L_0_i_beta(own.T);
   return polynom_p(u, i, j)*d2_G_alpha_wrt_c + (1 - polynom_p(u, i, j))*d2_G_beta_wrt_c;
}

__cuda_callable__
Real ZirconiumModel::deriv_of_g_w_resp_to_c_and_p(const VectorView& u, Index i, Index j) const
{
   Real c = conc_at(u, i, j);
   Real d_G_alpha_wrt_c = constants::G_Nb_alpha_0(own.T)
                            - constants::G_Zr_alpha_0(own.T)
                            + constants::R * own.T * (log(c) - log(1-c))
                            + (1 - 2*c) * constants::L_0_alpha;
   Real d_G_beta_wrt_c = constants::G_Nb_beta_0(own.T)
						   - constants::G_Zr_beta_0(own.T)
						   + constants::R * own.T * (log(c) - log(1-c))
						   + (1 - 2 * c) * constants::L_0_beta(own.T)
						   + (6*c - 6*pow(c, Real(2)) - 1) * constants::L_0_i_beta(own.T);

   return der_polynom_p(u, i, j)*(d_G_alpha_wrt_c - d_G_beta_wrt_c);
}

__cuda_callable__
Real ZirconiumModel::get_rhs_phase_at(const VectorView& u, Index i, Index j, Real F) const
{
    return Real(1.0)/own.alpha * (div_T0(u, i, j) + f_0(u, i , j) / own.ksi / own.ksi - own.par_b/own.ksi*grade_4_polynom(u, i, j)*F);
}

__cuda_callable__
Real ZirconiumModel::grad_p_1_forward(const VectorView& u, Index i, Index j) const
{
    return (phase_at(u, i+1, j) - phase_at(u, i, j))/grid.hx;
}

__cuda_callable__
Real ZirconiumModel::grad_p_2_forward(const VectorView& u, Index i, Index j) const
{
   return (phase_at(u, i, j+1) - phase_at(u, i, j))/grid.hy;
}

__cuda_callable__
Real ZirconiumModel::div_T0(const VectorView& u, Index i, Index j) const
{
    return (T0_1(grad_p_1_forward(u, i, j), grad_p_2_forward(u, i, j)) - T0_1(grad_p_1_forward(u, i-1, j), grad_p_2_forward(u, i-1, j)))/grid.hx
           + (T0_2(grad_p_1_forward(u, i, j), grad_p_2_forward(u, i, j)) - T0_2(grad_p_1_forward(u, i, j-1), grad_p_2_forward(u, i, j-1)))/grid.hy;
}

__cuda_callable__
Real ZirconiumModel::T0_1(const Real grad_p_1, const Real grad_p_2) const
{
    Real theta = atan2(-grad_p_2, -grad_p_1);
    const Anisotropy& a = own.anisotropy;
    return a.psi(theta) * a.psi(theta) * grad_p_1 - a.psi(theta) * a.der_psi(theta) * grad_p_2;
}

__cuda_callable__
Real ZirconiumModel::T0_2(const Real grad_p_1, const Real grad_p_2) const
{
    Real theta = atan2(-grad_p_2, -grad_p_1);
    const Anisotropy& a = own.anisotropy;
    return a.psi(theta) * a.psi(theta) * grad_p_2 + a.psi(theta) * a.der_psi(theta) * grad_p_1;
}

__cuda_callable__
Real ZirconiumModel::f_0(const VectorView& u, Index i, Index j) const
{
   return own.par_a*phase_at(u, i, j)*(1 - phase_at(u, i, j))*(phase_at(u, i, j) - Real(1.0)/Real(2.0));
}

__cuda_callable__
Real ZirconiumModel::grade_4_polynom(const VectorView& u, Index i, Index j) const
{
   return pow(phase_at(u, i, j), Real(2)) * pow(phase_at(u, i, j) - Real(1.0), Real(2));
}
