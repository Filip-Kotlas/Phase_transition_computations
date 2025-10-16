#pragma once
#include <cmath>
#include "types.hpp"

struct Domain {
    double x_left;
    double x_right;
    double y_left;
    double y_right;
};

enum class MODEL {
    MODEL_1 = 1,
    MODEL_2 = 2,
    MODEL_3 = 3,
    MODEL_4 = 4
};

enum class ICType {
        HyperbolicTangent,
        LinearByParts,
        ConstantCircle,
        ConstantHalves,
        Stripe,
        TwoBumps,
        Star,
        FourierX,
        FourierY,
        ThreeBumps,
        Count 
    };

namespace constants {
    inline constexpr Real R = 8.31446261815324;
    inline constexpr Real Cal_to_joule = 4.184;

    inline constexpr Real delta = 5e-9;
    inline constexpr Real delta_fixed = 5e-9;
    inline constexpr Real molar_volume = 1.4060e-5;
    inline constexpr Real b = 3.23e-10;
    inline constexpr Real sigma = 0.3;
    inline constexpr Real l = 50*delta;

    enum class Phase {
        alpha = 1,
        beta = 0
    };

    inline constexpr Real c_min = 0.0001;
    inline constexpr Real c_max = 1 - c_min;

    inline constexpr Real c_init_alpha = 0.007;
    inline constexpr Real c_init_beta = 0.025;

    inline constexpr Real L_0_alpha = 24411;
    __cuda_callable__ inline Real L_0_beta(Real T) {
        return 15911 + 3.35*T;
    }
    __cuda_callable__ inline Real L_0_i_beta(Real T){
        return 3919 - 1.091*T;
    }

    __cuda_callable__ inline Real G_Nb_alpha_0(Real T){
        return 1480.647
               + 144.445475*T
               - 26.4711*T*log(T)
               + 2.03475e-4*T*T
               - 3.5012e-7*T*T*T
               + 93399/T;
    }
    __cuda_callable__ inline Real G_Zr_alpha_0(Real T){
        return -7827.595
        + 125.64905*T
        - 24.1618*T*log(T)
        - 4.37791e-3*T*T
        + 34971/T;
    }
    __cuda_callable__ inline Real G_Nb_beta_0(Real T){
        return -8519.353
        + 142.045475*T
        - 26.4711*T*log(T)
        + 2.03475e-4*T*T
        - 3.5012e-7*T*T*T
        + 93399/T;
    }
    __cuda_callable__ inline Real G_Zr_beta_0(Real T){
        return -525.539
        + 124.9457*T
        - 25.607406*T*log(T)
        - 3.40084e-4*T*T
        - 9.729e-9*T*T*T
        + 25233/T
        - 7.6143e-11*T*T*T*T;
    }

    __cuda_callable__ inline Real G_m_alpha(Real c, Real T){
        return c*G_Nb_alpha_0(T)
        + (1-c)*G_Zr_alpha_0(T)
        + R*T*(c*log(c) + (1-c)*log(1-c))
        + c*(1-c)*L_0_alpha;
    }
    __cuda_callable__ inline Real G_m_beta(Real c, Real T){
        return c*G_Nb_beta_0(T)
        + (1-c)*G_Zr_beta_0(T)
        + R*T*(c*log(c) + (1-c)*log(1-c))
        + c*(1-c)*(L_0_beta(T) + L_0_i_beta(T)*(2*c-1));
    }

    __cuda_callable__ inline Real D_Nb_alpha(Real T){
        return 6.6e-10*exp(-31500*Cal_to_joule/(R*T));
    }
    __cuda_callable__ inline Real D_Nb_beta(Real T){
        return 9e-10*pow(T/1136.0, 18.1)*exp(-(25100+35.5*(T-1136))*Cal_to_joule/(R*T));
    }
    __cuda_callable__ inline Real D_eff(Real T){
        return pow(sqrt(D_Nb_alpha(T)) + sqrt(D_Nb_beta(T)), 2);
    }

    __cuda_callable__ inline Real M_Nb_alpha(Real T){
        return D_Nb_alpha(T) / R / T;
    }
    __cuda_callable__ inline Real M_Nb_beta(Real T){
        return D_Nb_beta(T) / R / T;
    }

    __cuda_callable__ inline Real M_phi(Real T){
        return 0.235 * D_eff(T) * delta_fixed / delta* molar_volume / (b*b*R*T);
    }
    __cuda_callable__ inline Real M_phi_tilde(Real T){
        return M_phi(T)*l*l/(M_Nb_alpha(T)*molar_volume);
    }

    __cuda_callable__ inline Real epsilon() {
        return sqrt(6*sigma*delta);
    }
    __cuda_callable__ inline Real epsilon_tilde(Real T){
        return epsilon()/(l*sqrt(R*T/molar_volume));
    }

    inline constexpr Real w = 3 * sigma / delta;
    __cuda_callable__ inline Real w_tilde(Real T){
        return w / (R*T);
    }

}