#ifndef ADAAI_DIFFERENTIATION_HPP
#define ADAAI_DIFFERENTIATION_HPP

#include <cmath>
#include <climits>

namespace ADAAI {
    enum class DiffMethod : int {
        stencil3,
        stencil3Extra,
        stencil5,
        stencil5Extra,
        FwdAAD
    };
    enum class WhichDerivative : int {
        x, // df/dx
        y, // df/dy
        xx, // d^2f/(dx)^2
        yy, // d^2f/(dy)^2
        xy // d^2f/dxdy
    };
    template<WhichDerivative W, DiffMethod M, typename Callable>
    double Differentiator(Callable const &F, double x, double y);

} // namespace ADAAI
#endif // ADAAI_DIFFERENTIATION_HPP