#include <cmath>

#include "ZirconiumModel.hpp"
#include "constants.hpp"

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
    return constants::G_m_alpha(conc_at(u, i, j), params.T) - constants::G_m_beta(conc_at(u, i, j), params.T);
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
    return params.par_b*params.par_d/Real(30.0)
            * conc_at(u, i, j)
            * (1 - conc_at(u, i, j))
            * constants::M_Nb_beta(params.T)
		    * pow(constants::M_Nb_alpha(params.T)/constants::M_Nb_beta(params.T), polynom_p(u, i, j))
		    * sec_deriv_of_g_w_resp_to_c(u, i, j);
}

__cuda_callable__
Real ZirconiumModel::get_phas_diff_coef(const VectorView& u, Index i, Index j) const
{
   return params.par_b*params.par_d/Real(30.0)
          * conc_at(u, i, j)
		  * (1 - conc_at(u, i, j))
          * constants::M_Nb_beta(params.T)
		  * pow(constants::M_Nb_alpha(params.T)/constants::M_Nb_beta(params.T), polynom_p(u, i, j))
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
   Real d2_G_alpha_wrt_c = constants::R * params.T / (c*(1-c))
   							 - 2 * constants::L_0_alpha;
   Real d2_G_beta_wrt_c = constants::R * params.T / (c*(1-c))
   							 - 2 * constants::L_0_beta(params.T)
							 + (6 - 12*c) * constants::L_0_i_beta(params.T);
   return polynom_p(u, i, j)*d2_G_alpha_wrt_c + (1 - polynom_p(u, i, j))*d2_G_beta_wrt_c;
}

__cuda_callable__
Real ZirconiumModel::deriv_of_g_w_resp_to_c_and_p(const VectorView& u, Index i, Index j) const
{
   Real c = conc_at(u, i, j);
   Real d_G_alpha_wrt_c = constants::G_Nb_alpha_0(params.T)
                            - constants::G_Zr_alpha_0(params.T)
                            + constants::R * params.T * (log(c) - log(1-c))
                            + (1 - 2*c) * constants::L_0_alpha;
   Real d_G_beta_wrt_c = constants::G_Nb_beta_0(params.T)
						   - constants::G_Zr_beta_0(params.T)
						   + constants::R * params.T * (log(c) - log(1-c))
						   + (1 - 2 * c) * constants::L_0_beta(params.T)
						   + (6*c - 6*pow(c, Real(2)) - 1) * constants::L_0_i_beta(params.T);

   return der_polynom_p(u, i, j)*(d_G_alpha_wrt_c - d_G_beta_wrt_c);
}
