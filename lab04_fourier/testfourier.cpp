#include <iostream>
#include "fourier.cpp"

int main() {
    int N = 5; // Количество членов в ряде Фурье

    for (int n = 0; n <= N; ++n) {
        double a_n = calculateFourierCoefficientA(n);
        double b_n = calculateFourierCoefficientB(n);
        std::cout << "a_" << n << " = " << a_n << ", b_" << n << " = " << b_n << std::endl;
    }

    return 0;
}