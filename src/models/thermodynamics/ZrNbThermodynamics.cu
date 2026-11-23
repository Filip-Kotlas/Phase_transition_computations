#include <cmath>

#include "models/thermodynamics/ZrNbThermodynamics.hpp"
#include "core/constants.hpp"

namespace {
// The functions of the temperature T [K], evaluated once, on the host, in double.
namespace from_T {
    using constants::R;
    using constants::Cal_to_joule;

    double L_0_beta(double T) {
        return 15911 + 3.35*T;
    }
    double L_0_i_beta(double T) {
        return 3919 - 1.091*T;
    }

    double G_Nb_alpha_0(double T) {
        return 1480.647
               + 144.445475*T
               - 26.4711*T*std::log(T)
               + 2.03475e-4*T*T
               - 3.5012e-7*T*T*T
               + 93399/T;
    }
    double G_Zr_alpha_0(double T) {
        return -7827.595
               + 125.64905*T
               - 24.1618*T*std::log(T)
               - 4.37791e-3*T*T
               + 34971/T;
    }
    double G_Nb_beta_0(double T) {
        return -8519.353
               + 142.045475*T
               - 26.4711*T*std::log(T)
               + 2.03475e-4*T*T
               - 3.5012e-7*T*T*T
               + 93399/T;
    }
    double G_Zr_beta_0(double T) {
        return -525.539
               + 124.9457*T
               - 25.607406*T*std::log(T)
               - 3.40084e-4*T*T
               - 9.729e-9*T*T*T
               + 25233/T
               - 7.6143e-11*T*T*T*T;
    }

    // the diffusivities of Nb in the two phases
    double D_Nb_alpha(double T) {
        return 6.6e-10*std::exp(-31500*Cal_to_joule/(R*T));
    }
    double D_Nb_beta(double T) {
        return 9e-10*std::pow(T/1136.0, 18.1)*std::exp(-(25100 + 35.5*(T - 1136))*Cal_to_joule/(R*T));
    }

    double M_Nb_alpha(double T) {
        return D_Nb_alpha(T) / R / T;
    }
    double M_Nb_beta(double T) {
        return D_Nb_beta(T) / R / T;
    }

    // Not used by the solver; kept for reference.
    // constexpr double delta = 5e-9;
    // constexpr double delta_fixed = 5e-9;
    // constexpr double molar_volume = 1.4060e-5;
    // constexpr double b = 3.23e-10;
    // constexpr double sigma = 0.3;
    // constexpr double l = 50*delta;
    // constexpr double w = 3*sigma/delta;
    // constexpr double c_min = 0.0001;
    // constexpr double c_max = 1 - c_min;
    // double D_eff(double T) {
    //     return std::pow(std::sqrt(D_Nb_alpha(T)) + std::sqrt(D_Nb_beta(T)), 2);
    // }
    // double M_phi(double T) {
    //     return 0.235*D_eff(T)*delta_fixed/delta*molar_volume/(b*b*R*T);
    // }
    // double M_phi_tilde(double T) {
    //     return M_phi(T)*l*l/(M_Nb_alpha(T)*molar_volume);
    // }
    // double epsilon() {
    //     return std::sqrt(6*sigma*delta);
    // }
    // double epsilon_tilde(double T) {
    //     return epsilon()/(l*std::sqrt(R*T/molar_volume));
    // }
    // double w_tilde(double T) {
    //     return w/(R*T);
    // }
}
}

ZrNbThermodynamics::ZrNbThermodynamics(double T)
    : RT(constants::R*T),
      G_Nb_alpha_0(from_T::G_Nb_alpha_0(T)),
      G_Zr_alpha_0(from_T::G_Zr_alpha_0(T)),
      G_Nb_beta_0(from_T::G_Nb_beta_0(T)),
      G_Zr_beta_0(from_T::G_Zr_beta_0(T)),
      L_0_beta(from_T::L_0_beta(T)),
      L_0_i_beta(from_T::L_0_i_beta(T)),
      M_Nb_alpha(from_T::M_Nb_alpha(T)),
      M_Nb_beta(from_T::M_Nb_beta(T))
{
}

__cuda_callable__
Real ZrNbThermodynamics::G_m_alpha(Real c) const
{
    return c*G_Nb_alpha_0
           + (1-c)*G_Zr_alpha_0
           + RT*(c*log(c) + (1-c)*log(1-c))
           + c*(1-c)*L_0_alpha;
}

__cuda_callable__
Real ZrNbThermodynamics::G_m_beta(Real c) const
{
    return c*G_Nb_beta_0
           + (1-c)*G_Zr_beta_0
           + RT*(c*log(c) + (1-c)*log(1-c))
           + c*(1-c)*(L_0_beta + L_0_i_beta*(2*c-1));
}
