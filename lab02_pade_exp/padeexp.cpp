#include "consts.hpp"
#include <cmath>

namespace ADAAI {
template <typename F> constexpr F PadeExp(F x) {
  static_assert(std::is_floating_point_v<F>);
  F numerator, denominator, pade_result;
  if (x < -c_Inf<F>) {
    return 0.0;
  }
  numerator = 1 + 2 * x / 5 + std::pow(x, 2) / 15 + std::pow(x, 3) / 180 +
              std::pow(x, 4) / 5040;
  denominator = 1 - 3 * x / 5 + std::pow(x, 2) / 6 - std::pow(x, 3) / 36 +
                std::pow(x, 4) / 336 - std::pow(x, 5) / 5040 +
                std::pow(x, 6) / 151200;
  pade_result = numerator / denominator;
  // (1+(2x)/5+x^{2}/15+x^{3}/180+x^{4}/5040)/(1-(3x)/5+x^{2}/6-x^{3}/36+x^{4}/336-x^{5}/5040+x^{6}/151200)
  // this is the closest to e^x pade approximant (built using wolfram)
  return pade_result;
}
} // namespace ADAAI
// вывести график экспоненты и паде из питона
// понять какой лучше брать чтобы была маленькая ошибка DONE
// выводить погрешность -- тесты DONE