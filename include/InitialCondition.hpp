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

        default:
            return static_cast<Real>(constants::Phase::beta);
            break;
        }
    }

    __cuda_callable__
    Real get_concentration(Index index) const
    {
        return 0.007 + (static_cast<Real>(constants::Phase::alpha) - get_phase(index)) * (0.025-0.007);
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