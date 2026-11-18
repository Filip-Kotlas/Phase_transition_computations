#include "InitialCondition.hpp"
#include <vector>
#include <iostream>
#include <TNL/Algorithms/parallelFor.h>
#include "ResultIO.hpp"

__cuda_callable__
Real InitialCondition::get_phase(Index ind) const
{
    Index i = ind % grid.sizeX;
    Index j = ind / grid.sizeX;
    Real x_coord = grid.x(i);
    Real y_coord = grid.y(j);
    Real x_rel = i * grid.hx; //relative to domain left
    Real y_rel = j * grid.hy; //relative to domain bottom
    Real x_central = x_coord - (grid.domain.x_right + grid.domain.x_left)/2;
    Real y_central = y_coord - (grid.domain.y_right + grid.domain.y_left)/2;

    Real dist_to_center = sqrt(pow(x_central, Real(2)) + pow(y_central, Real(2)));

    switch (init_condition)
    {
    case ICType::HyperbolicTangent:
        return Real(1.0)/2 * tanh(-3/ksi*(dist_to_center - init_cond_radius)) + Real(1.0)/2; //needs to be changed according to the correct phases
        break;

    case ICType::Circle:
        if( dist_to_center < init_cond_radius - Real(0.5)*slope_width )
        {
            return static_cast<Real>(constants::Phase::alpha);
        }
        else if( dist_to_center < init_cond_radius + Real(0.5)*slope_width )
        {
            return static_cast<Real>(constants::Phase::alpha) - (static_cast<Real>(constants::Phase::alpha) - static_cast<Real>(constants::Phase::beta))*(dist_to_center - (init_cond_radius - Real(0.5)*slope_width)) / (slope_width);
        }
        else
        {
            return static_cast<Real>(constants::Phase::beta);
        }
        break;

    case ICType::ConstantCircle:
        if( dist_to_center < init_cond_radius )
        {
            return static_cast<Real>(constants::Phase::alpha);
        }
        else
        {
            return static_cast<Real>(constants::Phase::beta);
        }
        break;

    case ICType::ConstantHalves:
        if( grid.domain.x_left + x_rel < (grid.domain.x_left + grid.domain.x_right) / 2 )
        {
            return static_cast<Real>(constants::Phase::alpha);
        }
        else
        {
            return static_cast<Real>(constants::Phase::beta);
        }
        break;

    case ICType::Stripe:
        if( x_rel < Real(0.2) )
        {
            return static_cast<Real>(constants::Phase::alpha);
        }
        else
        {
            return static_cast<Real>(constants::Phase::beta);
        }
        break;

    case ICType::TwoBumps: {
        Real y = 2*j*grid.hy;
        if( (y < 1 && i*grid.hx < y*y*(1-2*y+y*y)/Real(0.625)+Real(0.1)) ||
            (y >= 1 && i*grid.hx < (y-1)*(y-1)*(1-2*(y-1)+(y-1)*(y-1))/Real(0.625)+Real(0.1)))
        {
            return static_cast<Real>(constants::Phase::alpha);
        }
        else
        {
            return static_cast<Real>(constants::Phase::beta);
        }
        break;
    }

    case ICType::Star: {
        Real phi = atan2( y_central, x_central);
            
        if ( dist_to_center < Real(0.15) + Real(0.1) * sin(6 * phi) )
        {
            return static_cast<Real>(constants::Phase::alpha);
        }
        else
        {
            return static_cast<Real>(constants::Phase::beta);
        }
        break;
        }

    case ICType::ThreeBumps: {
        if( (pow(x_rel-Real(0.2), Real(2)) + pow(y_rel-Real(0.25), Real(2)) < Real(0.1)*Real(0.1)) ||
            (pow(x_rel-Real(0.2), Real(2)) + pow(y_rel-Real(0.5), Real(2)) < Real(0.1)*Real(0.1)) ||
            (pow(x_rel-Real(0.2), Real(2)) + pow(y_rel-Real(0.75), Real(2)) < Real(0.1)*Real(0.1)))
        {
            return static_cast<Real>(constants::Phase::alpha);
        }
        else
        {
            return static_cast<Real>(constants::Phase::beta);
        }
        break;
    }

    case ICType::PerpendicularStripes: {
        if( (grid.domain.x_left + x_rel < grid.domain.x_left + Real(0.1))
            || (grid.domain.y_left + y_rel > grid.domain.y_right - Real(0.1)) )
        {
            return static_cast<Real>(constants::Phase::alpha);
        }
        else
        {
            return static_cast<Real>(constants::Phase::beta);
        }
        break;
    }

    case ICType::Box:
        if( (grid.domain.x_left + x_rel < grid.domain.x_left + Real(0.1)
            || grid.domain.x_left + x_rel > grid.domain.x_right - Real(0.1))
            || (grid.domain.y_left + y_rel < grid.domain.y_left + Real(0.1)
            || grid.domain.y_left + y_rel > grid.domain.y_right - Real(0.1)) )
        {
            return static_cast<Real>(constants::Phase::alpha);
        }
        else
        {
            return static_cast<Real>(constants::Phase::beta);
        }
        break;

    case ICType::RandomBumps: {
        // Radius of the circles
        const Real cr = static_cast<Real>(0.03);

        // Fixed x coordinate of all circles
        const Real cx = static_cast<Real>(0.1);

        // y-interval for the placement (same range as originally: 0.1 .. 0.9)
        const Real y_min = static_cast<Real>(grid.domain.y_left + Real(0.03));
        const Real y_max = static_cast<Real>(grid.domain.y_right - Real(0.03));
        const Real y_range = y_max - y_min;


        // Reproducible seed (change it for a different layout or pass it from the host)
        const uint32_t seed = 1337u;

        bool inside_circle = false;

        Index number_of_circles = (grid.domain.y_right - grid.domain.y_left) * 10;

        for (Index n = 0; n < number_of_circles; ++n) {
            // t ∈ (n/20 .. (n+1)/20) with random "jitter"
            const Real t = (static_cast<Real>(n) + rnd01(static_cast<uint32_t>(n), seed)) / static_cast<Real>(number_of_circles);
            const Real cy = y_min + t * y_range;

            const Real dx = x_rel - cx;
            const Real dy = y_rel - cy;

            if (dx * dx + dy * dy < cr * cr) {
                inside_circle = true;
                break;
            }
        }

        if (inside_circle) {
            return static_cast<Real>(constants::Phase::alpha);
        } else {
            return static_cast<Real>(constants::Phase::beta);
        }
        break;
    }

    case ICType::WulffShape: {
        Real theta_origin = atan2(y_central, x_central);
        Real x_wulff = 0;
        Real y_wulff = 0;

        auto f = [&](Real wulff_theta) mutable {
            x_wulff = init_cond_radius*(anisotropy.psi(wulff_theta)*cos(wulff_theta) - anisotropy.der_psi(wulff_theta)*sin(wulff_theta));
            y_wulff = init_cond_radius*(anisotropy.psi(wulff_theta)*sin(wulff_theta) + anisotropy.der_psi(wulff_theta)*cos(wulff_theta));
            Real theta_computed = atan2(y_wulff, x_wulff);
            if (wulff_theta > constants::pi/2 && theta_computed < 0)
                theta_computed += 2*constants::pi;
            else if (wulff_theta < -constants::pi/2 && theta_computed > 0)
                theta_computed -= 2*constants::pi;
            return theta_computed - theta_origin;
        };

        auto compute_wulff_coords = [&]() -> void {
            Real wulff_left = -constants::pi;
            Real wulff_right = constants::pi;
            Real wulff_mid = (wulff_left + wulff_right)/2;

            // Declare f_left and f_right to check for root at the boundaries
            Real f_right = f(wulff_right);
            if (f_right == 0) {
                return;
            }
            Real f_left = f(wulff_left);
            if (f_left == 0) {
                return;
            }

            int iter_num = 0;

            // check if theta_origin is outside the range of function values at the boundaries and adjust accordingly
            if (f_left*f_right > 0) {
                theta_origin += (f_left > 0) ? 2*constants::pi : -2*constants::pi;
            }
            f_left = f(wulff_left);
            f_right = f(wulff_right);

            while(abs(f(wulff_mid)) > Real(1e-5)) {
                iter_num++;
                if (f_left*f(wulff_mid) < 0) {
                    wulff_right = wulff_mid;
                    wulff_mid = (wulff_left + wulff_right)/2;
                }
                else {
                    wulff_left =  wulff_mid;
                    wulff_mid = (wulff_left + wulff_right)/2;
                }
                if (iter_num > 100) {
                    printf("Wulff shape bisection did not converge! i: %d j: %d\n", i, j);
                    break;
                }
            }
        };
        
        compute_wulff_coords();
        Real radius_wulff = sqrt(x_wulff*x_wulff + y_wulff*y_wulff);
        if ( dist_to_center < radius_wulff - Real(0.5)*slope_width ) {
            return static_cast<Real>(constants::Phase::alpha);
        }
        else if ( dist_to_center < radius_wulff + Real(0.5)*slope_width ) {
            return static_cast<Real>(constants::Phase::alpha) - (static_cast<Real>(constants::Phase::alpha) - static_cast<Real>(constants::Phase::beta))*(dist_to_center - (radius_wulff - Real(0.5)*slope_width)) / (slope_width);
        }
        else {
            return static_cast<Real>(constants::Phase::beta);
        }
        break;
    }

    default:
        return static_cast<Real>(constants::Phase::beta);
        break;
    }
}

__cuda_callable__
Real InitialCondition::get_concentration(Index index) const {
    return constants::c_init_alpha
           + (static_cast<Real>(constants::Phase::alpha) - get_phase(index)) * (constants::c_init_beta - constants::c_init_alpha);
}

__cuda_callable__
uint32_t InitialCondition::wang_hash(uint32_t x) {
    x = (x ^ 61u) ^ (x >> 16);
    x *= 9u;
    x ^= (x >> 4);
    x *= 0x27d4eb2du;
    x ^= (x >> 15);
    return x;
}

__cuda_callable__
Real InitialCondition::rnd01(uint32_t idx, uint32_t seed) {
    uint32_t h = wang_hash(idx ^ seed);
    // 24-bit fraction -> (0,1)
    return (h & 0x00FFFFFFu) / static_cast<Real>(16777216.0);
}

void set_init_cond_manually(Vector& u, const InitialCondition& init_cond, const Grid& grid)
{
    const Index N = grid.num_points();

    u.forElements( 0, N,
        [ = ] __cuda_callable__( Index ind, Real& value )
        {
            value = init_cond.get_phase(ind);
        } );

    u.forElements( N, 2 * N,
        [ = ] __cuda_callable__( Index ind, Real& value )
        {
            value = init_cond.get_concentration(ind - N);
        } );
}

bool set_init_cond_from_file(Vector& u, const std::filesystem::path& filename, const Grid& grid, bool scaling)
{
    const Index N = grid.num_points();

    // Values from the file, the point (i, j) is on the position k = grid.index(i, j).
    std::vector<Real> file_phase;
    std::vector<Real> file_conc;
    if (!read_frame(filename, grid, file_phase, file_conc))
        return false;

    TNL::Containers::Vector<Real, TNL::Devices::Host, Index> host_u(2 * N);

    for (Index k = 0; k < N; ++k)
    {
        host_u[k] = static_cast<Real>(constants::Phase::beta); // phase initialization
        host_u[N + k] = constants::c_init_beta; // concentration initialization
    }

    const Index start_i = grid.sizeX / 4;
    const Index start_j = grid.sizeY / 4;

    for (Index count = 0; count < N; ++count) {
        if (scaling) {
            Index i = count % grid.sizeX;
            Index j = count / grid.sizeX;
            if (i % 2 == 0 && j % 2 == 0) {
                Index scale_index = grid.index(start_i + i / 2, start_j + j / 2);
                host_u[scale_index] = file_phase[count];
                host_u[N + scale_index] = file_conc[count];
            }
        }
        else {
            host_u[count] = file_phase[count];
            host_u[N + count] = file_conc[count];
        }
    }

    u = host_u;

    return true;
}

// Sets the boundary point edge of u from its neighbour (Neumann) or to dirichlet_value (Dirichlet).
__cuda_callable__ static
void set_boundary_point(VectorView& u, Index edge, Index neighbour, BCType bc, Real dirichlet_value)
{
    if (bc == BCType::Dirichlet)
        u[edge] = dirichlet_value;
    else if (bc == BCType::Neumann)
        u[edge] = u[neighbour];
}

void set_init_cond_boundary(Vector& u, const Grid& grid, const BoundaryConditions& bound_cond)
{
    const Index N = grid.num_points();
    const Real phase_dirichlet = static_cast<Real>(constants::Phase::beta);
    const Real conc_dirichlet = constants::c_init_beta;

    // Bottom and top edges first, then left and right edges, in the same order as the boundary
    // conditions in main.cu, so that the corners get the same values as there: when both
    // directions are Neumann, the corner becomes the value of the diagonal interior neighbour.
    auto bottom_top = [ = ] __cuda_callable__( Index i, VectorView& u )
    {
        const Index bottom = grid.index(i, 0);
        const Index top = grid.index(i, grid.sizeY - 1);
        const Index above_bottom = grid.index(i, 1);
        const Index below_top = grid.index(i, grid.sizeY - 2);
        set_boundary_point(u, bottom, above_bottom, bound_cond.phase_x, phase_dirichlet);
        set_boundary_point(u, top, below_top, bound_cond.phase_x, phase_dirichlet);
        set_boundary_point(u, N + bottom, N + above_bottom, bound_cond.conc_x, conc_dirichlet);
        set_boundary_point(u, N + top, N + below_top, bound_cond.conc_x, conc_dirichlet);
    };

    auto left_right = [ = ] __cuda_callable__( Index j, VectorView& u )
    {
        const Index left = grid.index(0, j);
        const Index right = grid.index(grid.sizeX - 1, j);
        const Index right_of_left = grid.index(1, j);
        const Index left_of_right = grid.index(grid.sizeX - 2, j);
        set_boundary_point(u, left, right_of_left, bound_cond.phase_y, phase_dirichlet);
        set_boundary_point(u, right, left_of_right, bound_cond.phase_y, phase_dirichlet);
        set_boundary_point(u, N + left, N + right_of_left, bound_cond.conc_y, conc_dirichlet);
        set_boundary_point(u, N + right, N + left_of_right, bound_cond.conc_y, conc_dirichlet);
    };

    VectorView view = u.getView();
    TNL::Algorithms::parallelFor< Device >(0, grid.sizeX, bottom_top, view);
    TNL::Algorithms::parallelFor< Device >(0, grid.sizeY, left_right, view);
}

bool set_initial_condition(Vector& u, const Parameters& parameters, Real ksi, const Anisotropy& anisotropy,
                           const std::filesystem::path& calculations_dir)
{
    const Grid grid = parameters.grid();

    if (parameters.init_cond_from_file){
        // The path in the config is relative to the calculations/ folder of this run.
        const std::filesystem::path init_cond_path = calculations_dir / parameters.init_cond_file_path;
        std::cout << "Setting initial condition from file: " << init_cond_path << std::endl;
        if (!set_init_cond_from_file(u, init_cond_path, grid, parameters.init_cond_scaling)){
            std::cerr << "Error: could not set the initial condition from the file " << init_cond_path
                      << " (init_cond_file_path = " << parameters.init_cond_file_path
                      << ", resolved relative to " << calculations_dir << "). Exiting." << std::endl;
            return false;
        }
    }
    else{
        std::cout << "Setting initial condition from code: " << ICType_to_string(parameters.init_condition) << std::endl;
        InitialCondition init_cond(parameters, ksi, anisotropy);
        set_init_cond_manually(u, init_cond, grid);
    }
    // The time stepping keeps the boundary values consistent only if they are consistent at the start.
    set_init_cond_boundary(u, grid, parameters.bound_cond);
    std::cout << "Initial condition set." << std::endl;
    return true;
}
