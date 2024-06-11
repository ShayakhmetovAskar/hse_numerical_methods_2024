#ifndef EVERHART_EVERHART_HPP
#define EVERHART_EVERHART_HPP
#include "consts.hpp"
#include <cmath>
#include <vector>

namespace ADAAI {

struct State {
  double x;
  double v;
};

double lagrange_basis(const std::array<double, ADAAI::N> &times, size_t i,
                      double t);

double lagrange_basis_derivative(const std::array<double, ADAAI::N> &times,
                                 size_t i, double t);

std::array<double, ADAAI::N> compute_coefficients(const State &state);

State interpolate(const std::array<State, ADAAI::N> &states,
                  const std::array<double, ADAAI::N> &times, double t);

State evolve_state(const State &state,
                   const std::array<double, ADAAI::N> &coeffs, double dt);

} // namespace ADAAI
#endif // EVERHART_EVERHART_HPP
