#include <iostream>
#include "fourier.cpp"

int main() {
    int N = 5; // Количество членов в ряде Фурье
    for (int n = 0; n <= N; ++n) { // to print coeffs
        if (n == 0) {
            double a_n = calculateFourierCoefficientA(n);
            std::cout << "a_" << n << " = " << a_n / 2 << ", b_" << n << " = " << 0 << '\n';
            n++;
        }
        double a_n = calculateFourierCoefficientA(n);
        double b_n = calculateFourierCoefficientB(n);
        std::cout << "a_" << n << " = " << a_n << ", b_" << n << " = " << b_n << '\n';
    }
    for (int n = 0; n < N; ++n) { //to print fourier series formula
        if (n == 0) {
            std::cout << calculateFourierCoefficientA(n) / 2;
            n++;
        }
        if(calculateFourierCoefficientA(n) < 0){
            std::cout << " " << calculateFourierCoefficientA(n) << " * cos(" << n << "x)";
        }
        else {
            std::cout << " + " << calculateFourierCoefficientA(n) << " * cos(" << n << "x)";
        }
        if (calculateFourierCoefficientB(n) < 0){
            std::cout << " " << calculateFourierCoefficientB(n) << " * sin(" << n << "x)";
        }
        else {
            std::cout << " + " << calculateFourierCoefficientB(n) << " * sin(" << n << "x)";
        }
    }
    return 0;
}