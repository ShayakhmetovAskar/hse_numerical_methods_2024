#include "consts.hpp"
#include "padeexp.cpp"
#include <cmath>
#include <iostream>

template <typename F> F checker() {
  F error;
  F max_error = 0;
  for (F x = ADAAI::c_Ln2<F> / (-2); x <= ADAAI::c_Ln2<F> / 2; x += 0.05) {
    error = std::abs(ADAAI::PadeExp(x) - std::exp(x));
    if (error >= max_error) {
      max_error = error;
    }
  }
  return max_error;
}

int main() {
  std::cout << checker<float>() << '\n';
  std::cout << checker<double>() << '\n';
  std::cout << checker<long double>() << '\n';
}
