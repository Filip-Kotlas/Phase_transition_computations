#pragma once
#include <cmath>

namespace constants {
    inline constexpr double R = 8.31446261815324;

    inline constexpr double L_0_alpha = 24411;
    inline double L_0_beta(int T) {
        return 15911 + 3.35*T;
    }
    inline double L_0_i_beta(int T){
        return 3919 - 1.091*T;
    }

    inline double G_Nb_alpha_0(int T){
        return 1480.647 + 144.445475*T - 26.4711*T*log(T) + 2.03475e-4*T*T - 3.5012e-7*T*T*T + 93399/T;
    }
    inline double G_Zr_alpha_0(int T){
        return -7827.595 + 125.64905*T - 24.1618*T*log(T) - 4.37791e-3*T*T + 34971/T;
    }
    inline double G_Nb_beta_0(int T){
        return -8519.353 + 142.045475*T - 26.4711*T*log(T) + 2.03475e-4*T*T - 3.5012e-7*T*T*T + 93399/T;
    }
    inline double G_Zr_beta_0(int T){
        return -525.539 + 124.9457*T - 25.607406*T*log(T) - 3.40084e-4*T*T - 9.729e-9*T*T*T + 25233/T - 7.6143e-11*T*T*T*T;
    }


    inline double D_Nb_alpha(int T){
        return 3;
        //return 6.6e-10*exp(-15851.4/T);
    }
    inline double D_Nb_beta(int T){
        return 1;
        //return 9e-9*pow(T/1136, 18.1)*exp(-(25100+35.5*(T-1136))/(1.98*T));
    }
    inline double D_eff(int T){
        return pow(sqrt(D_Nb_alpha(T)) + sqrt(D_Nb_beta(T)), 2);
    }

    inline double M_Nb_alpha(int T){
        return D_Nb_alpha(T) / R / T;
    }
    inline double M_Nb_beta(int T){
        return D_Nb_beta(T) / R / T;
    }

}
