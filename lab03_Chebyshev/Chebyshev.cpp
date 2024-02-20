#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <cmath>

namespace ADAAI {
    template<typename F>
    struct Chebyshev_series {
        F *coeffs; // coefficients
        int N; // order
        F lower_bound;
        F upper_bound;
        // x in [lower_bound, upper_bound]
        F *f; // Chebyshev function
    };

    template<typename F>
    F Chebyshev_count(Chebyshev_series<F> *cs, int N_input, F x) { // redefine a function from official documentation
        // to work with float as well.
        // used to calculate the value of Chebyshev polynomial in point x
        F T_n_minus_1 = 0.0;
        F T_n_minus_2 = 0.0;
        F m = (2.0 * x - cs->lower_bound - cs->upper_bound) / (cs->upper_bound - cs->lower_bound);
        F n = 2.0 * m;
        for (int i = N_input; i >= 1; i--) {
            F tmp = T_n_minus_1;
            T_n_minus_1 = n * T_n_minus_1 - T_n_minus_2 + cs->coeffs[i];
            T_n_minus_2 = tmp;
        }
        return m * T_n_minus_1 - T_n_minus_2 + 0.5 * cs->coeffs[0];
    }

    template<typename F>
    int Chebyshev_initialize(Chebyshev_series<F> *cs, gsl_function *func,
                             F lower_bound, F upper_bound) { // redefine a function from official documentation
        // to work with float as well
        cs->lower_bound = lower_bound;
        cs->upper_bound = upper_bound;

        for (int i = 0; i <= cs->N; i++) {
            F cheb_value = cos(M_PI * (i + 0.5) / (cs->N + 1));
            cs->f[i] = GSL_FN_EVAL(func, (cheb_value * 0.5 * (cs->upper_bound - cs->lower_bound) + 0.5 * (cs->upper_bound + cs->lower_bound)));
        }
        for (int j = 0; j <= cs->N; j++) {
            F coeff_sum = 0.0;
            for (int k = 0; k <= cs->N; k++) {
                coeff_sum += cs->f[k] * cos(M_PI * j * (k + 0.5) / (cs->N + 1));
            }
            cs->coeffs[j] = (2.0 / (cs->N + 1.0)) * coeff_sum;
        }
        return GSL_SUCCESS;
    }


    template<typename F>
    F Chebyshev(F x, int N) {
        static_assert(std::is_floating_point_v<F>);
        auto *cs = (Chebyshev_series<F> *) (malloc(sizeof(Chebyshev_series<F>)));
        cs->N = N;
        cs->coeffs = (F *) malloc((N + 1) * sizeof(F));
        cs->f = (F *) malloc((N + 1) * sizeof(F));
        gsl_function Chebyshev_approx; // gsl_function is declared in gsl_math.h
        Chebyshev_approx.function = [](double x, void *params) { return exp(x); }; //from official site
        Chebyshev_approx.params = nullptr; // from official site

        Chebyshev_initialize<F>(cs, &Chebyshev_approx, -1.0,
                                1.0); // initialize a Chebyshev function on the interval from -1.0 to 1.0

        F result = Chebyshev_count<F>(cs, N, x);
        free(cs);
        return result;
    }
}