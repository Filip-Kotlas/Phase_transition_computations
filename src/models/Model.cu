#include "models/Model.hpp"

Model::Model(const Parameters& parameters)
:  grid(parameters.grid()),
   bound_cond(parameters.bound_cond)
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

    if (bound_cond.phase_x == BCType::Dirichlet) {
        fu[bottom] = 0;
        fu[top] = 0;
    }
    else if (bound_cond.phase_x == BCType::Neumann) {
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

    if (bound_cond.phase_y == BCType::Dirichlet) {
        fu[left] = 0;
        fu[right] = 0;
    }
    else if (bound_cond.phase_y == BCType::Neumann) {
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

    if (bound_cond.conc_x == BCType::Dirichlet) {
        fu[bottom] = 0;
        fu[top] = 0;
    }
    else if (bound_cond.conc_x == BCType::Neumann) {
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

    if (bound_cond.conc_y == BCType::Dirichlet) {
        fu[left] = 0;
        fu[right] = 0;
    }
    else if (bound_cond.conc_y == BCType::Neumann) {
        fu[left] = fu[offset + grid.index(1, ind)];
        fu[right] = fu[offset + grid.index(grid.sizeX - 2, ind)];
    }
}
