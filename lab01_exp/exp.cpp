#include <climits>
#include <cmath>
#include <type_traits>
#include "exp.hpp"
#include "consts.hpp"


namespace ADAAI {
    template<typename F>
    constexpr F Exp(F x) {
        static_assert(std::is_floating_point_v<F>);
        F y = NAN;
        y = x / ln2<F>;
        F y0 = NAN;
        F y1 = std::modf(y, &y0); //y1 - fractional part, y0 - integral part
        if (y1 > 0.5) {
            y0++;
            y1--;
        } else if (y1 < -0.5) {
            y0--;
            y1++;
        }
        int n = (int) y0;
        if (n < INT_MIN)
            return 0.0;
        else if (n > INT_MAX)
            return inf<F>;
        else {
            F delta = 10.0 * eps<F>;
            F x1 = y1 * ln2<F>;
            F e_x1 = NAN;
            if (x1 > 0) {
                for (int i = 0; i <= n; i++) {
                    e_x1 += pow(x1, i) / tgamma(i + 1); // tgamma(i+1) - to calculate i!
                    if (delta > e_x1 * sqrt2<F>)
                        break;
                }
                F e_x = NAN;
                e_x = std::ldexp(e_x1, n); //e^x = e^x1 * 2^n
                return e_x;
            } else if (x1 < 0) {
                for (int i = 0; i <= n; i++) {
                    e_x1 += pow(x1, i) / tgamma(i + 1); // tgamma(i+1) - to calculate i!
                    if (delta > e_x1)
                        break;
                }
                F e_x = NAN;
                e_x = std::ldexp(e_x1, n); //e^x = e^x1 * 2^n
                return e_x;
            }
        }
    }
}