#include <cmath>
#include <algorithm>

#include "AllenCahnModel.hpp"

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
    if (params.force_term_type == FTType::Constant)
        return params.force_term_size;
    else if (params.force_term_type == FTType::Radial)
    {
        Real mid_x = (grid.domain.x_right - grid.domain.x_left)/2;
        Real mid_y = (grid.domain.y_right - grid.domain.y_left)/2;
        Real r = sqrt(pow(i*grid.hx - mid_x, Real(2)) + pow(j*grid.hy - mid_y, Real(2)));
        return -2/std::max(r, Real(0.1));
    }
    else
        return 0;
}
