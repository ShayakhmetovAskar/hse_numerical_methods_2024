#pragma once

#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <utility>
#include <vector>

constexpr double sigma_max = 0.29;

double f_sigma(double tau) {
  return tau < 0.25 ? 0.21 : tau < 0.5 ? 0.229 : tau < 0.75 ? 0.29 : 0.25;
}

double f_risk(double tau) {
  return tau <= 0.25   ? 0.16
         : tau <= 0.5  ? 0.15
         : tau <= 0.75 ? 0.14
         : tau <= 1.0  ? 0.13
                       : 0.12;
}

double int_r(double tau) {
  return tau <= 0.25   ? 0.16 * tau
         : tau <= 0.5  ? 0.04 + 0.15 * (tau - 0.25)
         : tau <= 0.75 ? 0.0775 + 0.14 * (tau - 0.5)
                       : 0.1125 + 0.13 * (tau - 0.75);
}

double int_sigma(double tau) {
  return tau <= 0.25   ? 0.21 * 0.21 * tau
         : tau <= 0.5  ? 0.011025 + 0.229 * 0.229 * (tau - 0.25)
         : tau <= 0.75 ? 0.02413525 + 0.29 * 0.29 * (tau - 0.5)
                       : 0.04516025 + 0.25 * 0.25 * (tau - 0.75);
}

namespace ADAAI::Integrator {
struct RHS {
  constexpr static int N = 0;
  virtual void operator()(double current_time, const double *current_state,
                          double *rhs) const = 0;
};
}

namespace ADAAI::Integrator {
template <typename RHS> struct Observer {
  virtual bool operator()(double current_time,
                          const double current_state[RHS::N]) const = 0;
};
}

struct a_rhs : Integrator::RHS {
  static constexpr int N = 502, K = 100;
  static constexpr double tau_max = 1.0, S_max = K * std::exp(5 * sigma_max);

  void operator()(double t, const double *u, double *rhs) const override {
    double s = f_sigma(t), s2 = s * s;

    for (int i = 1; i < N - 1; ++i) {
      double p = u[i - 1], c = u[i], n = u[i + 1];
      if (i == 1)
        p = 0;
      if (i == N - 2)
        n = S_max - K * std::exp(-int_r(t));

      rhs[i] = f_risk(t) * i * (n - p) / 2.0 +
               s2 * i * i * (n - 2 * c + p) / 2.0 - f_risk(t) * c;
    }
  }
};

struct a_func {
  static double get_consts(double *state, double S_tau) {
    int i = 0;
    while (i * a_rhs::S_max / a_rhs::N <= S_tau) {
      i++;
    }
    i--;
    S_tau -= i * a_rhs::S_max / a_rhs::N;
    return state[i] * (1.0 - S_tau) + state[i + 1] * S_tau;
  }
  static void init_start(double *state) {
    for (int i = 0; i < a_rhs::N; i++) {
      state[i] = std::max(i * a_rhs::S_max / a_rhs::N - a_rhs::K, 0.0);
    }
  }
};

struct a_obs : Integrator::Observer<a_rhs> {
  explicit a_obs() {}
  bool operator()(double current_time,
                  const double current_state[a_rhs::N]) const override {
    return current_time < 1.0;
  }
};

namespace ADAAI::Integrator::Stepper {
template <typename RHS> class T_Stepper {
  const RHS *m_rhs;
  constexpr static int N = RHS::N;
  explicit T_Stepper(const RHS *rhs) : m_rhs(rhs) {}
  virtual std::pair<double, double>
  operator()(double current_state[N], double next_state[N], double current_time,
             double suggested_d_time) const = 0;
};
}

namespace ADAAI::Integrator::Stepper {
template <typename RHS> class RKF_Stepper : public T_Stepper<RHS> {
public:
  explicit RKF_Stepper(const RHS *rhs) : T_Stepper<RHS>(rhs) {}

  std::pair<double, double>
  operator()(double current_state[RHS::N], double next_state[RHS::N],
             double current_time,
             double suggested_d_time = 0.01) const override {
    double h = suggested_d_time;
    std::vector<std::vector<double>> ks(6, std::vector<double>(RHS::N));
    for (int i = 0; i < 6; ++i) {
      double cur[RHS::N];
      std::memcpy(cur, current_state, RHS::N * sizeof(double));
      for (int j = 0; j < i; ++j) {
        double buf[RHS::N];
        std::memcpy(buf, ks[j].data(), RHS::N * sizeof(double));
        mul(buf, B_K_L[i + 1][j + 1], RHS::N);
        add(cur, buf, RHS::N);
      }
      (*this->m_rhs)(current_time + A_K[i + 1] * h, cur, ks[i].data());
      mul(ks[i].data(), h, RHS::N);
    }
    double cur[RHS::N] = {0};
    for (int i = 0; i < 6; ++i) {
      double buf[RHS::N];
      std::memcpy(buf, ks[i].data(), RHS::N * sizeof(double));
      mul(buf, CT_K[i + 1], RHS::N);
      add(cur, buf, RHS::N);
    }
    double TE = 0;
    for (int i = 0; i < RHS::N; ++i) {
      TE += cur[i] * cur[i];
    }
    double eps = 1e-9;
    double new_step = 0.9 * h * std::pow(eps / TE, 0.1);
    if (TE > eps) {
      return (*this)(current_state, next_state, current_time, new_step);
    }
    std::memcpy(next_state, current_state, RHS::N * sizeof(double));
    for (int i = 0; i < 6; ++i) {
      double buf[RHS::N];
      std::memcpy(buf, ks[i].data(), RHS::N * sizeof(double));
      mul(buf, CH_K[i + 1], RHS::N);
      add(next_state, buf, RHS::N);
    }
    return {current_time + h, new_step};
  }

private:
  static void mul(double *data, double c, std::size_t size) {
    for (std::size_t i = 0; i < size; ++i) {
      data[i] *= c;
    }
  }

  static void add(double *lhs, double *rhs, std::size_t size) {
    for (std::size_t i = 0; i < size; ++i) {
      lhs[i] += rhs[i];
    }
  }
  const std::vector<double> A_K = {0.0, 0, 0.5, 0.5, 1, 2.0 / 3, 0.2};
  const std::vector<double> C_K = {0.0, 1.0 / 6, 0, 2.0 / 3, 1.0 / 6};
  const std::vector<double> CH_K = {0.0,      1.0 / 24,  0,          0,
                                    5.0 / 48, 27.0 / 56, 125.0 / 336};
  const std::vector<double> CT_K = {0.0,      0.125,      0,           2.0 / 3,
                                    1.0 / 16, -27.0 / 56, -125.0 / 336};
  const std::vector<std::vector<double>> RKF_Stepper::B_K_L = {
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.5, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.25, 0.25, 0.0, 0.0, 0.0},
      {0.0, 0.0, -1.0, 2.0, 0.0, 0.0},
      {0.0, 7.0 / 27.0, 10.0 / 27.0, 0.0, 1.0 / 27.0, 0.0},
      {0.0, 28.0 / 625.0, -0.2, 546.0 / 625.0, 54.0 / 625.0, -378.0 / 625.0},
  };
};
}