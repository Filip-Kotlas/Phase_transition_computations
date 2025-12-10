#include "InitialCondition.hpp"

__cuda_callable__
Real InitialCondition::get_phase(Index ind) const
{
    Index i = ind % param.sizeX;
    Index j = ind / param.sizeX;
    Real x_coord = param.domain.x_left + i * hx;
    Real y_coord = param.domain.y_left + j * hy;
    Real x_rel = i * hx; //relative to domain left
    Real y_rel = j * hy; //relative to domain bottom

    Real dist_to_center = sqrt(pow(i*hx - (param.domain.x_right - param.domain.x_left)/2, 2)
                       + pow(j*hy - (param.domain.y_right - param.domain.y_left)/2, 2));

    switch (param.init_condition)
    {
    case ICType::HyperbolicTangent:
        return 1.0/2 * tanh(-3/param.ksi*(dist_to_center - init_cond_radius)) + 1.0/2; //needs to be changed according to the correct phases
        break;

    case ICType::LinearByParts:
        if( dist_to_center < init_cond_radius - 0.5*slope_width )
        {
            return static_cast<Real>(constants::Phase::alpha);
        }
        else if( dist_to_center < init_cond_radius + 0.5*slope_width )
        {
            return static_cast<Real>(constants::Phase::alpha) - (static_cast<Real>(constants::Phase::alpha) - static_cast<Real>(constants::Phase::beta))*(dist_to_center - (init_cond_radius - 0.5*slope_width)) / (slope_width);
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
        if( x_rel < (param.domain.x_left + param.domain.x_right) / 2 )
        {
            return static_cast<Real>(constants::Phase::alpha);
        }
        else
        {
            return static_cast<Real>(constants::Phase::beta);
        }
        break;

    case ICType::Stripe:
        if( x_rel < 0.2 )
        {
            return static_cast<Real>(constants::Phase::alpha);
        }
        else
        {
            return static_cast<Real>(constants::Phase::beta);
        }
        break;

    case ICType::TwoBumps: {
        Real y = 2*j*hy;
        if( (y < 1 && i*hx < y*y*(1-2*y+y*y)/0.625+0.1) ||
            (y >= 1 && i*hx < (y-1)*(y-1)*(1-2*(y-1)+(y-1)*(y-1))/0.625+0.1))
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
        Real phi = atan2( (y_rel - (param.domain.y_right-param.domain.y_left)/2), (x_rel - (param.domain.x_right - param.domain.x_left)/2));
            
        if ( dist_to_center < 0.15 + 0.1 * sin(6 * phi) )
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
        if( (pow(x_rel-0.2, 2) + pow(y_rel-0.25, 2) < 0.1*0.1) ||
            (pow(x_rel-0.2, 2) + pow(y_rel-0.5, 2) < 0.1*0.1) ||
            (pow(x_rel-0.2, 2) + pow(y_rel-0.75, 2) < 0.1*0.1))
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
        if( (x_rel < param.domain.x_left + 0.1) || (y_rel > param.domain.y_right - 0.1) )
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
        if( (x_rel < param.domain.x_left + 0.1 || x_rel > param.domain.x_right - 0.1) ||
            (y_rel < param.domain.y_left + 0.1 || y_rel > param.domain.y_right - 0.1) )
        {
            return static_cast<Real>(constants::Phase::alpha);
        }
        else
        {
            return static_cast<Real>(constants::Phase::beta);
        }
        break;

    case ICType::RandomBumps: {
        // Poloměr kruhů
        const Real cr = static_cast<Real>(0.03);

        // Fixní x souřadnice všech kruhů
        const Real cx = static_cast<Real>(0.1);

        // y-interval pro rozmístění (stejný rozsah jako původně: 0.1 .. 0.9)
        const Real y_min = static_cast<Real>(param.domain.y_left + 0.03);
        const Real y_max = static_cast<Real>(param.domain.y_right - 0.03);
        const Real y_range = y_max - y_min;


        // Reprodukovatelný seed (změň pro jiné rozložení nebo předej z hostu)
        const uint32_t seed = 1337u;

        bool inside_circle = false;

        Index number_of_circles = (param.domain.y_right - param.domain.y_left) * 10;

        for (Index n = 0; n < number_of_circles; ++n) {
            // t ∈ (n/20 .. (n+1)/20) s náhodným „jitterem“
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
        Real x = i*hx - (param.domain.x_right - param.domain.x_left)/2;
        Real y = j*hy - (param.domain.y_right - param.domain.y_left)/2;
        Real theta_origin = atan2(y, x);
        Real x_wulff = 0;
        Real y_wulff = 0;

        auto f = [&](Real wulff_theta) mutable {
            x_wulff = init_cond_radius*(psi(wulff_theta)*cos(wulff_theta) - der_psi(wulff_theta)*sin(wulff_theta));
            y_wulff = init_cond_radius*(psi(wulff_theta)*sin(wulff_theta) + der_psi(wulff_theta)*cos(wulff_theta));
            Real theta_computed = atan2(y_wulff, x_wulff);
            if (wulff_theta > M_PI/2 && theta_computed < 0)
                theta_computed += 2*M_PI;
            else if (wulff_theta < -M_PI/2 && theta_computed > 0)
                theta_computed -= 2*M_PI;
            return theta_computed - theta_origin;
        };

        auto compute_wulff_coords = [&]() -> void {
            Real wulff_left = -M_PI;
            Real wulff_right = M_PI;
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
                theta_origin += (f_left > 0) ? 2*M_PI : -2*M_PI;
            }
            f_left = f(wulff_left);
            f_right = f(wulff_right);

            while(abs(f(wulff_mid)) > 1e-5) {
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
        if ( dist_to_center < radius_wulff ) {
            return static_cast<Real>(constants::Phase::alpha);
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
    // 24bit frakce -> (0,1)
    return (h & 0x00FFFFFFu) / static_cast<Real>(16777216.0);
}

__cuda_callable__
Real InitialCondition::psi(const Real theta) const
{
    return 1 + param.A*cos(param.m*(theta - param.theta_0));
}

__cuda_callable__
Real InitialCondition::der_psi(const Real theta) const
{
    return -param.A*param.m*sin(param.m*(theta - param.theta_0));
}