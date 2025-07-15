#pragma once
#include <cmath>

namespace constants {
    inline constexpr double R = 8.31446261815324;
    inline constexpr double Cal_to_joule = 4.184;

    inline constexpr double delta = 5e-9;
    inline constexpr double molar_volume = 1.4060e-5;
    inline constexpr double b = 3.23e-10;
    inline constexpr double sigma = 0.3;
    inline constexpr double l = 50*delta;

    inline constexpr double p_alpha = 1;
    inline constexpr double p_beta = 0;

    inline constexpr double c_init_alpha = 0.007;
    inline constexpr double c_init_beta = 0.025;

    inline constexpr double L_0_alpha = 24411;
    inline double L_0_beta(float T) {
        return 15911 + 3.35*T;
    }
    inline double L_0_i_beta(float T){
        return 3919 - 1.091*T;
    }

    inline double G_Nb_alpha_0(float T){
        return 1480.647
               + 144.445475*T
               - 26.4711*T*log(T)
               + 2.03475e-4*T*T
               - 3.5012e-7*T*T*T
               + 93399/T;
    }
    inline double G_Zr_alpha_0(float T){
        return -7827.595
        + 125.64905*T
        - 24.1618*T*log(T)
        - 4.37791e-3*T*T
        + 34971/T;
    }
    inline double G_Nb_beta_0(float T){
        return -8519.353
        + 142.045475*T
        - 26.4711*T*log(T)
        + 2.03475e-4*T*T
        - 3.5012e-7*T*T*T
        + 93399/T;
    }
    inline double G_Zr_beta_0(float T){
        return -525.539
        + 124.9457*T
        - 25.607406*T*log(T)
        - 3.40084e-4*T*T
        - 9.729e-9*T*T*T
        + 25233/T
        - 7.6143e-11*T*T*T*T;
    }

    inline double G_m_alpha(double c, float T){
        return c*G_Nb_alpha_0(T)
        + (1-c)*G_Zr_alpha_0(T)
        + R*T*(c*log(c) + (1-c)*log(1-c))
        + c*(1-c)*L_0_alpha;
    }
    inline double G_m_beta(double c, float T){
        return c*G_Nb_beta_0(T)
        + (1-c)*G_Zr_beta_0(T)
        + R*T*(c*log(c) + (1-c)*log(1-c))
        + c*(1-c)*(L_0_beta(T) + L_0_i_beta(T)*(2*c-1));
    }

    inline double D_Nb_alpha(float T){
        return 6.6e-10*exp(-31500*Cal_to_joule/(R*T));
    }
    inline double D_Nb_beta(float T){
        return 9e-10*pow(T/1136.0, 18.1)*exp(-(25100+35.5*(T-1136))*Cal_to_joule/(R*T));
    }
    inline double D_eff(float T){
        return pow(sqrt(D_Nb_alpha(T)) + sqrt(D_Nb_beta(T)), 2);
    }

    inline double M_Nb_alpha(float T){
        return D_Nb_alpha(T) / R / T;
    }
    inline double M_Nb_beta(float T){
        return D_Nb_beta(T) / R / T;
    }

    inline double M_phi(float T){
        return 0.235*D_eff(T)*molar_volume/(b*b*R*T);
    }
    inline double M_phi_tilde(float T){
        return M_phi(T)*l*l/(M_Nb_alpha(T)*molar_volume);
    }

    inline double epsilon = sqrt(6*sigma*delta);
    inline double epsilon_tilde(float T){
        return epsilon/(l*sqrt(R*T/molar_volume));
    }

    inline double w = 3 * sigma / delta;
    inline double w_tilde(float T){
        return w / (R*T);
    }

}
