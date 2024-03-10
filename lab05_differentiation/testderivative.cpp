#include <cmath>
#include <iostream>
#include "differentiation.cpp"

double func(double x, double y) {
//    return std::sin(x * y) + std::exp(3 * x) * std::cos(y / x);
//    return 2 * x * std::sin(y) + x;
//    return 2 * std::sin(y) + x * y;
    return x * y * 2 + x * x;
}

int main() {
//    double x = 2.0763142;
//    double y = -0.07765;
    double x = 2.0;
    double y = 3.0;
    ADAAI::AAD22 z(0, x);
    ADAAI::AAD22 w(1, y);
//    ADAAI::AAD22 v = ADAAI::AAD22(2.0) * w.sin(w) + z * w;
//    ADAAI::AAD22 v = ADAAI::AAD22(2.0) * z * w.sin(w) + z;
    ADAAI::AAD22 v = ADAAI::AAD22(2.0) * z * w + z * z;
    long double result;
    result = ADAAI::Differentiator<ADAAI::WhichDerivative::x, ADAAI::DiffMethod::stencil3>(&func, x, y);
    std::cout << "df/dx using stencil3: " << result << '\n';
    result = ADAAI::Differentiator<ADAAI::WhichDerivative::x, ADAAI::DiffMethod::stencil3Extra>(&func, x, y);
    std::cout << "df/dx using stencil3Extra: " << result << '\n';
    result = ADAAI::Differentiator<ADAAI::WhichDerivative::x, ADAAI::DiffMethod::stencil5>(&func, x, y);
    std::cout << "df/dx using stencil5: " << result << '\n';
    result = ADAAI::Differentiator<ADAAI::WhichDerivative::x, ADAAI::DiffMethod::stencil5Extra>(&func, x, y);
    std::cout << "df/dx using stencil5Extra: " << result << '\n';
    std::cout << "df/dx using FwdAAD: " << v.m_d1[0] << '\n';
    std::cout << '\n';

    result = ADAAI::Differentiator<ADAAI::WhichDerivative::y, ADAAI::DiffMethod::stencil3>(&func, x, y);
    std::cout << "df/dy using stencil3: " << result << '\n';
    result = ADAAI::Differentiator<ADAAI::WhichDerivative::y, ADAAI::DiffMethod::stencil3Extra>(&func, x, y);
    std::cout << "df/dy using stencil3Extra: " << result << '\n';
    result = ADAAI::Differentiator<ADAAI::WhichDerivative::y, ADAAI::DiffMethod::stencil5>(&func, x, y);
    std::cout << "df/dy using stencil5: " << result << '\n';
    result = ADAAI::Differentiator<ADAAI::WhichDerivative::y, ADAAI::DiffMethod::stencil5Extra>(&func, x, y);
    std::cout << "df/dy using stencil5Extra: " << result << '\n';
    std::cout << "df/dy using FwdAAD: " << v.m_d1[1] << '\n';
    std::cout << '\n';

    result = ADAAI::Differentiator<ADAAI::WhichDerivative::xx, ADAAI::DiffMethod::stencil3>(&func, x, y);
    std::cout << "d^2f/(dx)^2 using stencil3: " << result << '\n';
    result = ADAAI::Differentiator<ADAAI::WhichDerivative::xx, ADAAI::DiffMethod::stencil3Extra>(&func, x, y);
    std::cout << "d^2f/(dx)^2 using stencil3Extra: " << result << '\n';
    result = ADAAI::Differentiator<ADAAI::WhichDerivative::xx, ADAAI::DiffMethod::stencil5>(&func, x, y);
    std::cout << "d^2f/(dx)^2 using stencil5: " << result << '\n';
    result = ADAAI::Differentiator<ADAAI::WhichDerivative::xx, ADAAI::DiffMethod::stencil5Extra>(&func, x, y);
    std::cout << "d^2f/(dx)^2 using stencil5Extra: " << result << '\n';
    std::cout << "d^2f/(dx)^2 using FwdAAD: " << v.m_d2[0] << '\n';
    std::cout << '\n';

    result = ADAAI::Differentiator<ADAAI::WhichDerivative::yy, ADAAI::DiffMethod::stencil3>(&func, x, y);
    std::cout << "d^2f/(dy)^2 using stencil3: " << result << '\n';
    result = ADAAI::Differentiator<ADAAI::WhichDerivative::yy, ADAAI::DiffMethod::stencil3Extra>(&func, x, y);
    std::cout << "d^2f/(dy)^2 using stencil3Extra: " << result << '\n';
    result = ADAAI::Differentiator<ADAAI::WhichDerivative::yy, ADAAI::DiffMethod::stencil5>(&func, x, y);
    std::cout << "d^2f/(dy)^2 using stencil5: " << result << '\n';
    result = ADAAI::Differentiator<ADAAI::WhichDerivative::yy, ADAAI::DiffMethod::stencil5Extra>(&func, x, y);
    std::cout << "d^2f/(dy)^2 using stencil5Extra: " << result << '\n';
    std::cout << "d^2f/(dy)^2 using FwdAAD: " << v.m_d2[1] << '\n';
    std::cout << '\n';

    result = ADAAI::Differentiator<ADAAI::WhichDerivative::xy, ADAAI::DiffMethod::stencil3>(&func, x, y);
    std::cout << "d^2f/dxdy using stencil3: " << result << '\n';
    result = ADAAI::Differentiator<ADAAI::WhichDerivative::xy, ADAAI::DiffMethod::stencil3Extra>(&func, x, y);
    std::cout << "d^2f/dxdy using stencil3Extra: " << result << '\n';
    result = ADAAI::Differentiator<ADAAI::WhichDerivative::xy, ADAAI::DiffMethod::stencil5>(&func, x, y);
    std::cout << "d^2f/dxdy using stencil5: " << result << '\n';
    result = ADAAI::Differentiator<ADAAI::WhichDerivative::xy, ADAAI::DiffMethod::stencil5Extra>(&func, x, y);
    std::cout << "d^2f/dxdy using stencil5Extra: " << result << '\n';
    std::cout << "d^2f/dxdy using FwdAAD: " << v.m_d2[2] << '\n';


    return 0;

}