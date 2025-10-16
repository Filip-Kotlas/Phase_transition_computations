#pragma once
#include "types.hpp"
#include "constants.hpp"

class InitialCondition {
public:

    InitialCondition() = default;

    InitialCondition(ICType type, Domain domain, Index sizeX, Index sizeY, Real ksi)
    : type(type), domain(domain), sizeX(sizeX), sizeY(sizeY), ksi(ksi) {
        r = (domain.x_right - domain.x_left)/6.0;
        r1 = r - 0.5*ksi;
        r2 = r1 + ksi;
        hx = (domain.x_right - domain.x_left) / (sizeX - 1);
        hy = (domain.y_right - domain.y_left) / (sizeY - 1);
    }

    __cuda_callable__
    Real get_phase(Index ind) const
    {
        Index i = ind % sizeX;
        Index j = ind / sizeX;
        Real radius = sqrt(pow(i*hx - (domain.x_right - domain.x_left)/2, 2) + pow(j*hy - (domain.y_right-domain.y_left)/2, 2));

        switch (type)
        {
        case ICType::HyperbolicTangent:
            return 1.0/2 * tanh(-3/ksi*(radius - r1)) + 1.0/2; //needs to be changed according to the correct phases
            break;

        case ICType::LinearByParts:
            if( radius < r1 )
            {
                return static_cast<Real>(constants::Phase::alpha);
            }
            else if( radius < r2 )
            {
                return static_cast<Real>(constants::Phase::alpha) - (static_cast<Real>(constants::Phase::alpha) - static_cast<Real>(constants::Phase::beta))*(radius - r1) / (r2 - r1);
            }
            else
            {
                return static_cast<Real>(constants::Phase::beta);
            }
            break;

        case ICType::ConstantCircle:
            if( radius < r1 )
            {
                return static_cast<Real>(constants::Phase::alpha);
            }
            else
            {
                return static_cast<Real>(constants::Phase::beta);
            }
            break;

        case ICType::ConstantHalves:
            if( i < sizeX/2 )
            {
                return static_cast<Real>(constants::Phase::alpha);
            }
            else
            {
                return static_cast<Real>(constants::Phase::beta);
            }
            break;

        case ICType::Stripe:
            if( i*hx < 0.2 )
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
            Real phi = atan( (j*hy - (domain.y_right-domain.y_left)/2) / (i*hx - (domain.x_right - domain.x_left)/2));
            if ((i*hx - (domain.x_right - domain.x_left)/2) < 0 )
                phi = phi + M_PI;
            else if ( (i*hx - (domain.x_right - domain.x_left)/2) > 0 && (j*hy - (domain.y_right-domain.y_left)/2) < 0)
                phi = phi + 2 * M_PI;
            
            if ( radius < 0.15 + 0.1 * sin(6 * phi) )
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
            Real x = i*hx;
            Real y = j*hy;
            if( (pow(x-0.2, 2) + pow(y-0.25, 2) < 0.1*0.1) ||
                (pow(x-0.2, 2) + pow(y-0.5, 2) < 0.1*0.1) ||
                (pow(x-0.2, 2) + pow(y-0.75, 2) < 0.1*0.1))
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
            if( (i*hx < domain.x_left + 0.1) || (j*hy > domain.y_right - 0.1) )
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
            if( (i*hx < domain.x_left + 0.1 || i*hx > domain.x_right - 0.1) ||
                (j*hy < domain.y_left + 0.1 || j*hy > domain.y_right - 0.1) )
            {
                return static_cast<Real>(constants::Phase::alpha);
            }
            else
            {
                return static_cast<Real>(constants::Phase::beta);
            }
            break;

        case ICType::RandomBumps: {
            const Real x = i * hx;
            const Real y = j * hy;

            // Poloměr kruhů
            const Real cr = static_cast<Real>(0.03);

            // Fixní x souřadnice všech kruhů
            const Real cx = static_cast<Real>(0.1);

            // y-interval pro rozmístění (stejný rozsah jako původně: 0.1 .. 0.9)
            const Real y_min = static_cast<Real>(domain.y_left + 0.03);
            const Real y_max = static_cast<Real>(domain.y_right - 0.03);
            const Real y_range = y_max - y_min;


            // Reprodukovatelný seed (změň pro jiné rozložení nebo předej z hostu)
            const uint32_t seed = 1337u;

            bool inside_circle = false;

            Index number_of_circles = (domain.y_right - domain.y_left) * 10;

            // 20 kruhů, stratifikované vzorkování: rozsekáme y-interval na 20 binů a v každém náhodně posuneme střed
            for (Index n = 0; n < number_of_circles; ++n) {
                // t ∈ (n/20 .. (n+1)/20) s náhodným „jitterem“
                const Real t = (static_cast<Real>(n) + rnd01(static_cast<uint32_t>(n), seed)) / static_cast<Real>(number_of_circles);
                const Real cy = y_min + t * y_range;

                const Real dx = x - cx;
                const Real dy = y - cy;

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

        default:
            return static_cast<Real>(constants::Phase::beta);
            break;
        }
    }

    __cuda_callable__
    Real get_concentration(Index index) const {
        return 0.007 + (static_cast<Real>(constants::Phase::alpha) - get_phase(index)) * (0.025-0.007);
    }

    // Pomocné funkce: jednoduchý hash a uniformní [0,1)
    __cuda_callable__ static
    uint32_t wang_hash(uint32_t x) {
        x = (x ^ 61u) ^ (x >> 16);
        x *= 9u;
        x ^= (x >> 4);
        x *= 0x27d4eb2du;
        x ^= (x >> 15);
        return x;
    }

    __cuda_callable__ static
    Real rnd01(uint32_t idx, uint32_t seed) {
        uint32_t h = wang_hash(idx ^ seed);
        // 24bit frakce -> (0,1)
        return (h & 0x00FFFFFFu) / static_cast<Real>(16777216.0);
    }

private:
    ICType type;
    Domain domain;
    Index sizeX;
    Index sizeY;
    Real ksi;

    Real r;
    Real r1;
    Real r2;
    Real hx;
    Real hy;
};