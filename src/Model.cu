#include <cmath>

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
Real Model::get_rhs_phase_at(const VectorView& u, Index i, Index j, Real F) const
{
    return Real(1.0)/params.alpha * (div_T0(u, i, j) + f_0(u, i , j) / params.ksi / params.ksi - params.par_b/params.ksi*grade_4_polynom(u, i, j)*F);
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
Real Model::f_0(const VectorView& u, Index i, Index j) const
{
   return params.par_a*phase_at(u, i, j)*(1 - phase_at(u, i, j))*(phase_at(u, i, j) - Real(1.0)/Real(2.0));
}

__cuda_callable__
Real Model::grade_4_polynom(const VectorView& u, Index i, Index j) const
{
   return pow(phase_at(u, i, j), Real(2)) * pow(phase_at(u, i, j) - Real(1.0), Real(2));
}