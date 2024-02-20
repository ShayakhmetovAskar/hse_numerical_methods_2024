#include <gsl/gsl_chebyshev.h>
#include <cmath>

namespace ADAAI {
    template<typename F>
    F Chebyshev(F x, size_t order) {
        static_assert(std::is_floating_point_v<F>);
        gsl_cheb_series *cs = gsl_cheb_alloc(order); // allocate memory

        gsl_function Chebyshev_approx;
        Chebyshev_approx.function = [](double x, void *params) { return exp(x); };
        Chebyshev_approx.params = nullptr;

        gsl_cheb_init(cs, &Chebyshev_approx, -1.0, 1.0); // initialize a Chebyshev function on the interval from -1.0 to 1.0

        F result = gsl_cheb_eval_n(cs, order, x); // info from official site:
        // This function evaluates the Chebyshev series cs at a given point x, to (at most) the given order.

        gsl_cheb_free(cs); // free memory

        return result;
    }
}