#ifndef SATELLITE_SATELLITE_HPP
#define SATELLITE_SATELLITE_HPP
#include <iostream>

namespace ADAAI {
class EverhartIntegrator {
public:
  EverhartIntegrator(double step) {}
  void integrate(std::vector<double> &position, std::vector<double> &velocity,
                 int steps, const std::string &filename);

private:
  double dt;
  std::vector<double> alpha;
  std::vector<double> coeff_position;
  std::vector<double> coeff_velocity;

  void calculateCoefficients() {}
  std::vector<double> computeAcceleration(const std::vector<double> &position);
  void everhartInterpolation(std::vector<double> &position,
                             std::vector<double> &velocity);
};
} // namespace ADAAI

#endif // SATELLITE_SATELLITE_HPP
