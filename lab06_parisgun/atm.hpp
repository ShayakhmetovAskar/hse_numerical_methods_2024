#ifndef ADAAI_ATM_HPP
#define ADAAI_ATM_HPP
#include <vector>

namespace ADAAI {
class SAM { // SAM = Standard Atmosphere Model
public:
  std::pair<double, double> operator()(double h);
  static double get_pressure(double h);
  static double get_density(double h);
  static double get_A(double h);
};
} // namespace ADAAI

#endif // ADAAI_ATM_HPP
