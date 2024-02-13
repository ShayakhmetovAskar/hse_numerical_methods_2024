#include <cmath>
#include <iostream>

#include "exp.cpp"
#include "consts.hpp"

template<typename F>
F checker() {
    F error;
    F max_error = 0;
    for (F x = ADAAI::c_Ln2<F> / (-2); x <= ADAAI::c_Ln2<F> / 2; x += 0.05) {
        if (x <= 0) {
            error = std::abs(ADAAI::Exp(x) - std::exp(x));
        } else {
            error = std::abs(ADAAI::Exp(x) / std::exp(x) - 1.0);
        }
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