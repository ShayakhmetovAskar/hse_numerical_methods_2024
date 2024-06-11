#include "fourier.cpp"
#include <iostream>

int main() {
  int N = 5; // Количество членов в ряде Фурье
  double x[N + 1];
  double result = calculateFourierCoefficientA(0) / 2;
  for (int n = 0; n <= N; ++n) { // to print coeffs
    if (n == 0) {
      double a_n = calculateFourierCoefficientA(n);
      std::cout << "a_" << n << " = " << a_n / 2 << ", b_" << n << " = " << 0
                << '\n';
      n++;
    }
    double a_n = calculateFourierCoefficientA(n);
    double b_n = calculateFourierCoefficientB(n);
    std::cout << "a_" << n << " = " << a_n << ", b_" << n << " = " << b_n
              << '\n';
  }
  for (int n = 0; n <= N; n++) {
    if (n == 0) {
      x[0] = 0;
      std::cout << result << ' ';
      std::cout << std::exp(x[0]) << '\n';
      n++;
    }
    x[n] = M_PI * n / (N + 1);
    result += calculateFourierCoefficientA(n) * std::cos(n * x[n]) +
              calculateFourierCoefficientB(n) * std::sin(n * x[n]);
    std::cout << result << ' ';
    std::cout << std::exp(x[n]) << '\n';
  }
  for (int n = 0; n < N; ++n) { // to print fourier series formula
    if (n == 0) {
      std::cout << calculateFourierCoefficientA(n) / 2;
      n++;
    }
    if (calculateFourierCoefficientA(n) < 0) {
      std::cout << " " << calculateFourierCoefficientA(n) << " * cos(" << n
                << "x)";
    } else {
      std::cout << " + " << calculateFourierCoefficientA(n) << " * cos(" << n
                << "x)";
    }
    if (calculateFourierCoefficientB(n) < 0) {
      std::cout << " " << calculateFourierCoefficientB(n) << " * sin(" << n
                << "x)";
    } else {
      std::cout << " + " << calculateFourierCoefficientB(n) << " * sin(" << n
                << "x)";
    }
  }
  return 0;
}