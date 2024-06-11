#ifndef ADAAI_PG_HPP
#define ADAAI_PG_HPP
#include <tuple>

namespace ADAAI {
class ParisGun { // SAM = Standard Atmosphere Model
public:
  static const int n = 13;
  static constexpr inline double speed[n] = { // speed in machs
      0.4, 0.6, 0.8, 0.86, 0.91, 0.95, 1.0, 1.2, 1.4, 1.61, 1.8, 2.1, 2.2};
  static constexpr inline double c_d[n] = {0.086, 0.097, 0.1,  0.11,  0.125,
                                           0.215, 0.37,  0.32, 0.305, 0.285,
                                           0.27,  0.235, 0.22};
  static double interpolate(double v_M);
  static int projectile_equations(double t, const double u[], double dudt[],
                                  void *params);
  static void simulate_projectile_motion();
};
} // namespace ADAAI

#endif // ADAAI_PG_HPP
