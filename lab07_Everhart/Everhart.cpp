#include "Everhart.hpp"
#include "consts.hpp"
#include <array>
#include <cmath>
#include <iostream>
#include <limits>

namespace ADAAI {

double ADAAI::lagrange_basis(const std::array<double, 4> &times, size_t i,
                             double t) {
  double result = 1.0;
  for (size_t j = 0; j < 4; ++j) {
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

double ADAAI::lagrange_basis_derivative(const std::array<double, 4> &times,
                                        size_t i, double t) {
  double result = 0.0;
  for (size_t j = 0; j < 4; ++j) {
    if (j != i) {
      double term = 1.0 / (times[i] - times[j]);
      if (std::fabs(times[i] - times[j]) <
          std::numeric_limits<double>::epsilon()) {
        std::cerr << "Denominator too small in lagrange_basis_derivative\n";
        return 0.0;
      }
      for (size_t l = 0; l < 4; ++l) {
        if (l != i && l != j) {
          term *= (t - times[l]) / (times[i] - times[l]);
          if (std::fabs(times[i] - times[l]) <
              std::numeric_limits<double>::epsilon()) {
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

std::array<double, 4> ADAAI::compute_coefficients(const State &state) {
  std::array<double, 4> coeffs;
  coeffs[0] = state.x;
  coeffs[1] = state.v;
  coeffs[2] = -k * state.x;
  coeffs[3] = -k * state.v;

  for (size_t i = 2; i < 4; ++i) {
    coeffs[i] /= static_cast<double>(i * (i - 1));
  }

  return coeffs;
}

State ADAAI::interpolate(const std::array<State, 4> &states,
                         const std::array<double, 4> &times, double t) {
  double p = 0.0;
  double q = 0.0;
  for (size_t i = 0; i < 4; ++i) {
    double basis = ADAAI::lagrange_basis(times, i, t);
    double basis_derivative = ADAAI::lagrange_basis_derivative(times, i, t);
    p += states[i].x * basis;
    q += states[i].v * basis_derivative;
  }
  return {p, q};
}

State ADAAI::evolve_state(const State &state,
                          const std::array<double, 4> &coeffs, double dt) {
  State new_state = state;
  double dt_pow = 1.0;
  for (size_t i = 1; i < 4; ++i) {
    dt_pow *= dt;
    new_state.x += coeffs[i] * dt_pow;
    new_state.v += (i < 3 ? coeffs[i + 1] : 0) * dt_pow;
  }
  return new_state;
}
} // namespace ADAAI

int main() {
  std::array<ADAAI::State, 4> states = {
      ADAAI::State{1.0, 0.0}, ADAAI::State{0.995, -0.1},
      ADAAI::State{0.98, -0.2}, ADAAI::State{0.955, -0.3}};
  std::array<double, 4> times = {0.0, 1.0, 2.0, 3.0};

  double t = 4.0;

  for (int step = 0; step < 100; ++step) {
    ADAAI::State current_state = states[3];
    auto coeffs = ADAAI::compute_coefficients(current_state);
    ADAAI::State new_state = ADAAI::evolve_state(current_state, coeffs, 1.0);

    for (size_t i = 0; i < 3; ++i) {
      states[i] = states[i + 1];
      times[i] = times[i + 1];
    }
    states[3] = new_state;
    times[3] = t;
    t += 1.0;

    std::cout << "t = " << t << ", x = " << new_state.x
              << ", v = " << new_state.v << "\n";
  }

  return 0;
}