#include <cmath>
#include <iostream>
#include "Chebyshev.cpp"
#include "consts.hpp"


template<typename F>
F checker() {
    F error;
    F max_error = 0;
    for (F x = -1; x <= 1; x += 0.05) {
        int N = 20; // seems like order = 20 is the best, because it gives the smallest error.
        error = std::abs(ADAAI::Chebyshev(x, N) - std::exp(x));
        if (error >= max_error) {
            max_error = error;
        }
    }
    return max_error;
}

int main() {
    std::cout << checker<float>() << '\n'; // for order = 20: max_error is approximately equal to 2*eps
//    std::cout << ADAAI::c_Eps<float> << '\n';
    std::cout << checker<double>() << '\n'; // for order = 20: max_error is approximately equal to 10*eps
//    std::cout << ADAAI::c_Eps<double> << '\n';
    std::cout << checker<long double>() << '\n'; // for order = 20: max_error is approximately equal to 10*eps
//    std::cout << ADAAI::c_Eps<long double> << '\n';
}