#include <cmath>
#include <iostream>

#include "exp.cpp"
#include "consts.hpp"

template<typename F>
bool checker() {
    F error;
    F max_error = 0;
    for (F x = -50; x <= 50; x++) {
        if (x <= 0) {
            error = std::abs(ADAAI::Exp(x) - std::exp(x));
        } else {
            error = std::abs(ADAAI::Exp(x) / std::exp(x) - 1.0);
        }
        if (error >= max_error) {
            max_error = error;
        }
    }
    return error <= ADAAI::c_Eps<F> * 10;
}

int main() {
    std::cout << checker<double>() << '\n';
    std::cout << checker<float>() << '\n';
    std::cout << checker<long double>() << '\n';
}
