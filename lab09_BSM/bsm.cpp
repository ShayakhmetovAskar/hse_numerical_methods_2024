#include "bsm.hpp"
#include <cmath>
#include <cstdint>
#include <utility>

namespace ADAAI::Integration::PDE_BSM::Implicit {
template <typename RHS>
class ImplicitStepper : public Integration::Integrator::Stepper::T_Stepper<RHS> {
  constexpr static int N = RHS::N;

  double fill_A(uint64_t i, double tau, [[maybe_unused]] double dTau) const {
    return i * AUX_FUNC::f_sigma(tau) * i * AUX_FUNC::f_sigma(tau) / 2 - i * AUX_FUNC::f_risk(tau);
  }

  double fill_B(uint64_t i, double tau, double d_tau) const {
    return -AUX_FUNC::f_risk(tau) - i * AUX_FUNC::f_sigma(tau) * i * AUX_FUNC::f_sigma(tau) - 1 / d_tau;
  }

  double fill_D(uint64_t i, double tau, double d_tau) const {
    return i * AUX_FUNC::f_sigma(tau) * i * AUX_FUNC::f_sigma(tau) / 2 + i * AUX_FUNC::f_risk(tau);
  }

  void fill_matrix(double matrix[N][N], double tau, double d_tau) const {
    for (int i = 1; i <= N - 2; ++i) {
      matrix[i][i] = fill_B(i, tau, d_tau);
      matrix[i + 1][i] = fill_A(i + 1, tau, d_tau);
      matrix[i][i + 1] = fill_D(i, tau, d_tau);
    }
    matrix[N - 2][N - 2] = fill_B(N - 2, tau, d_tau);
  }

  double init_start_state() const { return 0.0; }

  double get_last_state(double tau) const {
    return a_rhs::S_max - a_rhs::K * std::exp(-AUX_FUNC::int_r(tau));
  }

  void fill_F(double F_i[N], double matrix[N][N], double current_state[RHS::N], double tau, double d_tau) const {
    for (uint64_t i = 1; i <= N - 2; ++i) {
      F_i[i] = -current_state[i] / d_tau;
    }

    F_i[1] -= fill_A(1, tau, d_tau) * get_zero_state();
    F_i[N - 2] -= fill_D(N - 2, tau, d_tau) * get_last_state(tau);

    for (int i = 1; i <= N - 3; ++i) {
      double coef = matrix[i + 1][i] / matrix[i][i];
      matrix[i + 1][i] = 0.0;
      matrix[i + 1][i + 1] -= coef * matrix[i][i + 1];
      F_i[i + 1] -= coef * F_i[i];
    }

    for (int i = N - 2; i >= 2; --i) {
      double coef = matrix[i - 1][i] / matrix[i][i];
      matrix[i - 1][i] = 0.0;
      F_i[i - 1] -= coef * F_i[i];
    }
  }

  void find_new_state(double F_i[N], double matrix[N][N], double next_state[RHS::N], double tau) const {
    for (int i = 1; i <= N - 2; ++i) {
      next_state[i] = F_i[i] / matrix[i][i];
    }

    next_state[0] = get_zero_state();
    next_state[N - 1] = get_last_state(tau);
  }

  explicit ImplicitStepper(RHS *rhs) : Integrator::Stepper::T_Stepper<RHS>(rhs) {}

  std::pair<double, double> operator()(double cur_state[RHS::N], double next_state[RHS::N], double cur_time, double m_d_time) const override {
    static double matrix[N][N];
    static double F_i[N];

    fill_matrix(matrix, cur_time, m_d_time);
    fill_F(F_i, matrix, cur_state, cur_time, m_d_time);
    find_new_state(F_i, matrix, next_state, cur_time);

    return {cur_time + m_d_time, m_d_time};
  }
};
}

namespace ADAAI::Integration::PDE_BSM::Numerical {
enum class e_sol { EXPLICIT, IMPLICIT };

double solve(double S_m_tau, double m_tau, e_sol sol) {
  double d_tau = m_tau / 1000;
  double state[a_rhs::N], end_state[a_rhs::N];
  ADAAI::PDE_BSM::init_start(state);

  auto rhs = a_rhs();
  auto observer = a_obs();

  if (sol == e_sol::EXPLICIT) {
    auto stepper = Integrator::Stepper::RKF_Stepper(&rhs);
    auto integrator = Integrator::ODE_Integrator<a_rhs, Integrator::Stepper::RKF_Stepper<a_rhs>>(&stepper, &observer);
    integrator(state, end_state, 0.0, m_tau, d_tau);
  } else {
    auto stepper = Implicit::ImplicitStepper(&rhs);
    auto integrator = Integrator::ODE_Integrator<a_rhs, Implicit::ImplicitStepper<a_rhs>>(&stepper, &observer);
    integrator(state, end_state, 0.0, m_tau, d_tau);
  }
  return a_func::get_consts(end_state, S_m_tau);
}
}