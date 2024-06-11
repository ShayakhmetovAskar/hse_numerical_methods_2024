#include <iostream>
#include <array>
#include <cmath>
#include <limits>
#include "Everhart.hpp"
#include "consts.hpp"


namespace ADAAI {

    double ADAAI::lagrange_basis(const std::array<double, N> &times, size_t i, double t) {
        double result = 1.0;
        for (size_t j = 0; j < N; ++j) {
            if (j != i) {
                double denom = times[i] - times[j];
                if (std::fabs(denom) < std::numeric_limits<double>::epsilon()) {
                    std::cerr << "Denominator too small in lagrange_basis\n";
                    return 0.0;
                }
                result *= (t - times[j]) / denom;
            }
        }
        return result;
    }

    double ADAAI::lagrange_basis_derivative(const std::array<double, N> &times, size_t i, double t) {
        double result = 0.0;
        for (size_t j = 0; j < N; ++j) {
            if (j != i) {
                double term = 1.0 / (times[i] - times[j]);
                if (std::fabs(times[i] - times[j]) < std::numeric_limits<double>::epsilon()) {
                    std::cerr << "Denominator too small in lagrange_basis_derivative\n";
                    return 0.0;
                }
                for (size_t l = 0; l < N; ++l) {
                    if (l != i && l != j) {
                        term *= (t - times[l]) / (times[i] - times[l]);
                        if (std::fabs(times[i] - times[l]) < std::numeric_limits<double>::epsilon()) {
                            std::cerr << "Denominator too small in lagrange_basis_derivative\n";
                            return 0.0;
                        }
                    }
                }
                result += term;
            }
        }
        return result;
    }

    std::array<double, N> ADAAI::compute_coefficients(const State &state) {
        std::array<double, N> coeffs;
        coeffs[0] = state.x;
        coeffs[1] = state.v;
        coeffs[2] = -k * state.x;
        coeffs[3] = -k * state.v;

        for (size_t i = 2; i < N; ++i) {
            coeffs[i] /= static_cast<double>(i * (i - 1));
        }

        return coeffs;
    }

    State ADAAI::interpolate(const std::array<State, N> &states, const std::array<double, N> &times, double t) {
        double p = 0.0;
        double q = 0.0;
        for (size_t i = 0; i < N; ++i) {
            double basis = lagrange_basis(times, i, t);
            double basis_derivative = lagrange_basis_derivative(times, i, t);
            p += states[i].x * basis;
            q += states[i].v * basis_derivative;
        }
        return {p, q};
    }

    State ADAAI::evolve_state(const State &state, const std::array<double, N> &coeffs, double dt) {
        State new_state = state;
        double dt_pow = 1.0;
        for (size_t i = 1; i < N; ++i) {
            dt_pow *= dt;
            new_state.x += coeffs[i] * dt_pow;
            new_state.v += (i < N - 1 ? coeffs[i + 1] : 0) * dt_pow;
        }
        return new_state;
    }
}
int main() {
    std::array<ADAAI::State, ADAAI::N> states = {
            ADAAI::State{1.0, 0.0},
            ADAAI::State{0.995, -0.1},
            ADAAI::State{0.98, -0.2},
            ADAAI::State{0.955, -0.3}
    };
    std::array<double, ADAAI::N> times = {0.0, ADAAI::dt, 2 * ADAAI::dt, 3 * ADAAI::dt};

    double t = 4 * ADAAI::dt;

    for (int step = 0; step < 100; ++step) {
        ADAAI::State current_state = states[ADAAI::N - 1];
        auto coeffs = ADAAI::compute_coefficients(current_state);
        ADAAI::State new_state = evolve_state(current_state, coeffs, ADAAI::dt);

        for (size_t i = 0; i < ADAAI::N - 1; ++i) {
            states[i] = states[i + 1];
            times[i] = times[i + 1];
        }
        states[ADAAI::N - 1] = new_state;
        times[ADAAI::N - 1] = t;
        t += ADAAI::dt;

        std::cout << "t = " << t << ", x = " << new_state.x << ", v = " << new_state.v << "\n";
    }

    return 0;
}