#include <type_traits>
#include <cmath>
#include <climits>

template<typename F>
constexpr inline F ln2;
template<>
constexpr inline float ln2<float> = 0.6931472f;
template<>
constexpr inline double ln2<double> = 0.693147180559945;
template<>
constexpr inline long double ln2<long double> = 0.6931471805599453L;

template<typename F>
constexpr inline F eps;
template<>
constexpr inline float eps<float> = std::numeric_limits<float>::epsilon();
template<>
constexpr inline double eps<double> = std::numeric_limits<double>::epsilon();
template<>
constexpr inline long double eps<long double> = std::numeric_limits<long double>::epsilon();

template<typename F>
constexpr inline F sqrt2;
template<>
constexpr inline float sqrt2<float> = (float) M_SQRT2;
template<>
constexpr inline float sqrt2<double> = (double) M_SQRT2;
template<>
constexpr inline float sqrt2<long double> = (long double) M_SQRT2;

template<typename F>
constexpr F Exp(F x) {
    static_assert(std::is_floating_point_v<F>);
    F y;
    y = x / ln2<F>;
    F y0;
    F y1 = std::modf(y, &y0); //y1 - fractional part, y0 - integral part
    if (y1 > 0.5) {
        y0 += 1;
        y1 = y1 - 1;
    } else if (y1 < -0.5) {
        y0 -= 1;
        y1 = y1 + 1;
    }
    int n = (int) y0;
    if (n < INT_MIN) {
        return 0.0;
    } else if (n > INT_MAX) {
        return std::numeric_limits<F>::infinity();
    } else {
        F delta = 10.0 * eps<F>;
        F x1 = y1 * ln2<F>;
        F e_x1;
        if (x1 > 0) {
            for (int i = 0; i <= n; i++) {
                e_x1 += pow(x1, i) / tgamma(i + 1); // tgamma(i+1) - to calculate i!
                if (delta > e_x1 * sqrt2<F>){
                    break;
                }
            }
            F e_x;
            e_x = std::ldexp(e_x1, n); //e^x = e^x1 * 2^n
            return e_x;
        }
        else if (x1 < 0){
            for (int i = 0; i <= n; i++) {
                e_x1 += pow(x1, i) / tgamma(i + 1); // tgamma(i+1) - to calculate i!
                if (delta > e_x1){
                    break;
                }
            }
            F e_x;
            e_x = std::ldexp(e_x1, n); //e^x = e^x1 * 2^n
            return e_x;
        }
    }
}

