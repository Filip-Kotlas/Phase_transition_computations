#pragma once

namespace constants {
    constexpr double R = 8.31446261815324;
    constexpr double L_0_alpha = 24411;

    double L_0_beta(int T) {
        return 15911 + 3.35*T;
    }
    double L_0_i_beta(int T){
        return 3919 - 1.091*T;
    }
    
}
