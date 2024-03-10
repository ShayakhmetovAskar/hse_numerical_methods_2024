#include <cmath>
#include <iostream>
#include "differentiation.cpp"

double func(double x, double y) {
    return std::sin(x * y) + std::exp(3 * x) * std::cos(y / x);
}

int main() {
    // works well for every x, y
    // example:
    double x = 2.0763142;
    double y = -0.07765;
    long double result;
    result = ADAAI::Differentiator<ADAAI::WhichDerivative::x, ADAAI::DiffMethod::stencil3>(&func, x, y);
    std::cout << "df/dx using stencil3: " << result << '\n';
    result = ADAAI::Differentiator<ADAAI::WhichDerivative::x, ADAAI::DiffMethod::stencil3Extra>(&func, x, y);
    std::cout << "df/dx using stencil3Extra: " << result << '\n';
    result = ADAAI::Differentiator<ADAAI::WhichDerivative::x, ADAAI::DiffMethod::stencil5>(&func, x, y);
    std::cout << "df/dx using stencil5: " << result << '\n';
    result = ADAAI::Differentiator<ADAAI::WhichDerivative::x, ADAAI::DiffMethod::stencil5Extra>(&func, x, y);
    std::cout << "df/dx using stencil5Extra: " << result << '\n';
//    result = ADAAI::Differentiator<ADAAI::WhichDerivative::x, ADAAI::DiffMethod::FwdAAD>(&func, x, y);
//    std::cout << "df/dx using FwdAAD" << result << '\n';
    std::cout << '\n';

    result = ADAAI::Differentiator<ADAAI::WhichDerivative::y, ADAAI::DiffMethod::stencil3>(&func, x, y);
    std::cout << "df/dy using stencil3: " << result << '\n';
    result = ADAAI::Differentiator<ADAAI::WhichDerivative::y, ADAAI::DiffMethod::stencil3Extra>(&func, x, y);
    std::cout << "df/dy using stencil3Extra: " << result << '\n';
    result = ADAAI::Differentiator<ADAAI::WhichDerivative::y, ADAAI::DiffMethod::stencil5>(&func, x, y);
    std::cout << "df/dy using stencil5: " << result << '\n';
    result = ADAAI::Differentiator<ADAAI::WhichDerivative::y, ADAAI::DiffMethod::stencil5Extra>(&func, x, y);
    std::cout << "df/dy using stencil5Extra: " << result << '\n';
//    result = ADAAI::Differentiator<ADAAI::WhichDerivative::y, ADAAI::DiffMethod::FwdAAD>(&func, x, y);
//    std::cout << "df/dy using FwdAAD" << result << '\n';
    std::cout << '\n';

    result = ADAAI::Differentiator<ADAAI::WhichDerivative::xx, ADAAI::DiffMethod::stencil3>(&func, x, y);
    std::cout << "d^2f/(dx)^2 using stencil3: " << result << '\n';
    result = ADAAI::Differentiator<ADAAI::WhichDerivative::xx, ADAAI::DiffMethod::stencil3Extra>(&func, x, y);
    std::cout << "d^2f/(dx)^2 using stencil3Extra: " << result << '\n';
    result = ADAAI::Differentiator<ADAAI::WhichDerivative::xx, ADAAI::DiffMethod::stencil5>(&func, x, y);
    std::cout << "d^2f/(dx)^2 using stencil5: " << result << '\n';
    result = ADAAI::Differentiator<ADAAI::WhichDerivative::xx, ADAAI::DiffMethod::stencil5Extra>(&func, x, y);
    std::cout << "d^2f/(dx)^2 using stencil5Extra: " << result << '\n';
//    result = ADAAI::Differentiator<ADAAI::WhichDerivative::xx, ADAAI::DiffMethod::FwdAAD>(&func, x, y);
//    std::cout << "d^2f/(dx)^2 using FwdAAD" << result << '\n';
    std::cout << '\n';

    result = ADAAI::Differentiator<ADAAI::WhichDerivative::yy, ADAAI::DiffMethod::stencil3>(&func, x, y);
    std::cout << "d^2f/(dy)^2 using stencil3: " << result << '\n';
    result = ADAAI::Differentiator<ADAAI::WhichDerivative::yy, ADAAI::DiffMethod::stencil3Extra>(&func, x, y);
    std::cout << "d^2f/(dy)^2 using stencil3Extra: " << result << '\n';
    result = ADAAI::Differentiator<ADAAI::WhichDerivative::yy, ADAAI::DiffMethod::stencil5>(&func, x, y);
    std::cout << "d^2f/(dy)^2 using stencil5: " << result << '\n';
    result = ADAAI::Differentiator<ADAAI::WhichDerivative::yy, ADAAI::DiffMethod::stencil5Extra>(&func, x, y);
    std::cout << "d^2f/(dy)^2 using stencil5Extra: " << result << '\n';
//    result = ADAAI::Differentiator<ADAAI::WhichDerivative::yy, ADAAI::DiffMethod::FwdAAD>(&func, x, y);
//    std::cout << "d^2f/(dy)^2 using FwdAAD" << result << '\n';
    std::cout << '\n';

    result = ADAAI::Differentiator<ADAAI::WhichDerivative::xy, ADAAI::DiffMethod::stencil3>(&func, x, y);
    std::cout << "d^2f/dxdy using stencil3: " << result << '\n';
    result = ADAAI::Differentiator<ADAAI::WhichDerivative::xy, ADAAI::DiffMethod::stencil3Extra>(&func, x, y);
    std::cout << "d^2f/dxdy using stencil3Extra: " << result << '\n';
    result = ADAAI::Differentiator<ADAAI::WhichDerivative::xy, ADAAI::DiffMethod::stencil5>(&func, x, y);
    std::cout << "d^2f/dxdy using stencil5: " << result << '\n';
    result = ADAAI::Differentiator<ADAAI::WhichDerivative::xy, ADAAI::DiffMethod::stencil5Extra>(&func, x, y);
    std::cout << "d^2f/dxdy using stencil5Extra: " << result << '\n';
//    result = ADAAI::Differentiator<ADAAI::WhichDerivative::xy, ADAAI::DiffMethod::FwdAAD>(&func, x, y);
//    std::cout << "d^2f/dxdy using FwdAAD" << result << '\n';


    return 0;

}