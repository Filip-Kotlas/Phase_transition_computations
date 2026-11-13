#include "Model.hpp"

Model::Model(const Parameters& parameters)
:  grid(parameters.grid()),
   params(parameters.model_params())
{}

Index Model::getDegreesOfFreedom() const
{
    return 2 * grid.num_points();
}

__cuda_callable__
void Model::set_rhs_at(const VectorView& u, VectorView& fu, Index i, Index j) const
{
    const Index offset = grid.num_points();
    fu[grid.index(i, j)] = get_rhs_phase_at(u, i, j);
    fu[offset + grid.index(i, j)] = get_rhs_concentration_at(u, i, j);
}

__cuda_callable__
void Model::apply_boundary_condition_xdir(Index ind, VectorView fu) const
{
    apply_phase_boundary_condition_xdir(ind, fu);
    apply_concentration_boundary_condition_xdir(ind, fu);
}

__cuda_callable__
void Model::apply_boundary_condition_ydir(Index ind, VectorView fu) const
{
    apply_phase_boundary_condition_ydir(ind, fu);
    apply_concentration_boundary_condition_ydir(ind, fu);
}

__cuda_callable__
void Model::apply_phase_boundary_condition_xdir(Index ind, VectorView fu) const
{
    // Bottom and top edges, ind = i
    const Index bottom = grid.index(ind, 0);
    const Index top = grid.index(ind, grid.sizeY - 1);

    if (params.bc_phase_x == BCType::Dirichlet) {
        fu[bottom] = 0;
        fu[top] = 0;
    }
    else if (params.bc_phase_x == BCType::Neumann) {
        fu[bottom] = fu[grid.index(ind, 1)];
        fu[top] = fu[grid.index(ind, grid.sizeY - 2)];
    }
}

__cuda_callable__
void Model::apply_phase_boundary_condition_ydir(Index ind, VectorView fu) const
{
    // Left and right edges, ind = j
    const Index left = grid.index(0, ind);
    const Index right = grid.index(grid.sizeX - 1, ind);

    if (params.bc_phase_y == BCType::Dirichlet) {
        fu[left] = 0;
        fu[right] = 0;
    }
    else if (params.bc_phase_y == BCType::Neumann) {
        fu[left] = fu[grid.index(1, ind)];
        fu[right] = fu[grid.index(grid.sizeX - 2, ind)];
    }
}

__cuda_callable__
void Model::apply_concentration_boundary_condition_xdir(Index ind, VectorView fu) const
{
    // Bottom and top edges, ind = i
    const Index offset = grid.num_points();
    const Index bottom = offset + grid.index(ind, 0);
    const Index top = offset + grid.index(ind, grid.sizeY - 1);

    if (params.bc_conc_x == BCType::Dirichlet) {
        fu[bottom] = 0;
        fu[top] = 0;
    }
    else if (params.bc_conc_x == BCType::Neumann) {
        fu[bottom] = fu[offset + grid.index(ind, 1)];
        fu[top] = fu[offset + grid.index(ind, grid.sizeY - 2)];
    }
}

__cuda_callable__
void Model::apply_concentration_boundary_condition_ydir(Index ind, VectorView fu) const
{
    // Left and right edges, ind = j
    const Index offset = grid.num_points();
    const Index left = offset + grid.index(0, ind);
    const Index right = offset + grid.index(grid.sizeX - 1, ind);

    if (params.bc_conc_y == BCType::Dirichlet) {
        fu[left] = 0;
        fu[right] = 0;
    }
    else if (params.bc_conc_y == BCType::Neumann) {
        fu[left] = fu[offset + grid.index(1, ind)];
        fu[right] = fu[offset + grid.index(grid.sizeX - 2, ind)];
    }
}

__cuda_callable__
Real Model::get_rhs_phase_at(const VectorView& u, Index i, Index j) const
{
    return Real(1.0)/params.alpha * (div_T0(u, i, j) + f_0(u, i , j) / params.ksi / params.ksi - params.par_b/params.ksi*grade_4_polynom(u, i, j)*F(u, i, j));
}

__cuda_callable__
Real Model::get_rhs_concentration_at(const VectorView& u, Index i, Index j) const
{
    if (params.force_term_type == FTType::Zirconium)
        return div_D_grad_concentration(u, i, j) + div_D_grad_phase(u, i, j);
    else
        return 0;
}

__cuda_callable__
Real Model::grad_p_1_forward(const VectorView& u, Index i, Index j) const
{
    return (phase_at(u, i+1, j) - phase_at(u, i, j))/grid.hx;
}

__cuda_callable__
Real Model::grad_p_2_forward(const VectorView& u, Index i, Index j) const
{
   return (phase_at(u, i, j+1) - phase_at(u, i, j))/grid.hy;
}

__cuda_callable__
Real Model::div_T0(const VectorView& u, Index i, Index j) const
{
    return (T0_1(grad_p_1_forward(u, i, j), grad_p_2_forward(u, i, j)) - T0_1(grad_p_1_forward(u, i-1, j), grad_p_2_forward(u, i-1, j)))/grid.hx
           + (T0_2(grad_p_1_forward(u, i, j), grad_p_2_forward(u, i, j)) - T0_2(grad_p_1_forward(u, i, j-1), grad_p_2_forward(u, i, j-1)))/grid.hy;
}

__cuda_callable__
Real Model::T0_1(const Real grad_p_1, const Real grad_p_2) const
{
    Real theta = atan2(-grad_p_2, -grad_p_1);
    const Anisotropy& a = params.anisotropy;
    return a.psi(theta) * a.psi(theta) * grad_p_1 - a.psi(theta) * a.der_psi(theta) * grad_p_2;
}

__cuda_callable__
Real Model::T0_2(const Real grad_p_1, const Real grad_p_2) const
{
    Real theta = atan2(-grad_p_2, -grad_p_1);
    const Anisotropy& a = params.anisotropy;
    return a.psi(theta) * a.psi(theta) * grad_p_2 + a.psi(theta) * a.der_psi(theta) * grad_p_1;
}

__cuda_callable__
Real Model::div_D_grad_concentration(const VectorView& u, Index i, Index j) const
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
Real Model::div_D_grad_phase(const VectorView& u, Index i, Index j) const
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
Real Model::get_conc_diff_coef(const VectorView& u, Index i, Index j) const
{
    return params.par_b*params.par_d/Real(30.0)
            * conc_at(u, i, j)
            * (1 - conc_at(u, i, j))
            * constants::M_Nb_beta(params.T)
		    * pow(constants::M_Nb_alpha(params.T)/constants::M_Nb_beta(params.T), polynom_p(u, i, j))
		    * sec_deriv_of_g_w_resp_to_c(u, i, j);
}

__cuda_callable__
Real Model::get_phas_diff_coef(const VectorView& u, Index i, Index j) const
{
   return params.par_b*params.par_d/Real(30.0)
          * conc_at(u, i, j)
		  * (1 - conc_at(u, i, j))
          * constants::M_Nb_beta(params.T)
		  * pow(constants::M_Nb_alpha(params.T)/constants::M_Nb_beta(params.T), polynom_p(u, i, j))
		  * deriv_of_g_w_resp_to_c_and_p(u, i, j);
}

__cuda_callable__
Real Model::f_0(const VectorView& u, Index i, Index j) const
{
   return params.par_a*phase_at(u, i, j)*(1 - phase_at(u, i, j))*(phase_at(u, i, j) - Real(1.0)/Real(2.0));
}

__cuda_callable__
Real Model::F(const VectorView& u, Index i, Index j) const
{
    if (params.force_term_type == FTType::Constant)
        return params.force_term_size;
    else if (params.force_term_type == FTType::Radial)
    {
        Real mid_x = (grid.domain.x_right - grid.domain.x_left)/2;
        Real mid_y = (grid.domain.y_right - grid.domain.y_left)/2;
        Real r = sqrt(pow(i*grid.hx - mid_x, Real(2)) + pow(j*grid.hy - mid_y, Real(2)));
        return -2/std::max(r, Real(0.1));
    }
    else if (params.force_term_type == FTType::Zirconium)
        return constants::G_m_alpha(conc_at(u, i, j), params.T) - constants::G_m_beta(conc_at(u, i, j), params.T);
    else
        return 0;

}

__cuda_callable__
Real Model::grade_4_polynom(const VectorView& u, Index i, Index j) const
{
   return pow(phase_at(u, i, j), Real(2)) * pow(phase_at(u, i, j) - Real(1.0), Real(2));
}

__cuda_callable__
Real Model::polynom_p(const VectorView& u, Index i, Index j) const
{
   return 6*pow(phase_at(u, i, j), Real(5)) - 15*pow(phase_at(u, i, j), Real(4)) + 10*pow(phase_at(u, i, j), Real(3));
}

__cuda_callable__
Real Model::der_polynom_p(const VectorView& u, Index i, Index j) const
{
   return 30*(pow(phase_at(u, i, j), Real(4)) - 2*pow(phase_at(u, i, j), Real(3)) + pow(phase_at(u, i, j), Real(2)));
}

__cuda_callable__
Real Model::sec_deriv_of_g_w_resp_to_c(const VectorView& u, Index i, Index j) const
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
Real Model::deriv_of_g_w_resp_to_c_and_p(const VectorView& u, Index i, Index j) const
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