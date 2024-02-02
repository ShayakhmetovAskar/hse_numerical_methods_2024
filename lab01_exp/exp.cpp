#include <climits>
#include <cmath>
#include "exp.hpp"
#include "consts.hpp"


namespace ADAAI {
    template<typename F>
    constexpr F Exp(F x) {
        static_assert(std::is_floating_point_v<F>);
        F y0, y1;
        y1 = std::modf(x, &y0);
        if (y1 > 0.5) {
            y0++;
            y1--;
        } if (y1 < -0.5) {
            y0--;
            y1++;
        }
        if (y0 < INT_MIN)
            return 0.0;
        if (y0 > INT_MAX)
            return c_Inf<F>;
        F x0 = 0, x1;
        F y = x / c_Ln2<F>;
        x1 = std::modf(y, &x0);
        if (x1 > 0.5) {
            x0++;
            x1--;
        } if (x1 < -0.5) {
            x0--;
            x1++;
        }
        if (x0 < INT_MIN)
            return 0.0;
        if (x0 > INT_MAX)
            return c_Inf<F>;
        F x2 = x1 * c_Ln2<F>;
        bool is_neg = false;
        if (x2 < 0)
            is_neg = true;
        F e_x1 = 0.0;
        long long n = 1;
        F taylor_part = 1.0;
        F delta = 10.0 * c_Eps<F>;
        while (taylor_part * c_Sqrt2<F> >= delta) {
            e_x1 += taylor_part;
            taylor_part *= x2;
            taylor_part /= ++n;
        }
        F e_x = std::ldexp(e_x1, static_cast<int>(x0));
        return is_neg ? 1/e_x : e_x;
    }
}