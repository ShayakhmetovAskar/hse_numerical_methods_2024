#include "differentiation.hpp"
#include <cmath>
#include <functional>

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
double Differentiator(Callable const &F, double x, double y) {
    double h_x, h_y;
    int n = 2; // can be changed to any integer > 2 - is used to make mistake smaller
    if (abs(x) < 1) {
        h_x = 1e-4;
    } else {
        h_x = abs(x) * 1e-4;
    }
    if (abs(y) < 1) {
        h_y = 1e-4;
    } else {
        h_y = abs(y) * 1e-4;
    }
    if (M == DiffMethod::stencil3) {
        if constexpr (W == WhichDerivative::x) {
            double df_dx = (F(x + h_x, y) - F(x - h_x, y)) / (2 * h_x);
            return df_dx;
        } else if (W == WhichDerivative::y) {
            double df_dy = (F(x, y + h_y) - F(x, y - h_y)) / (2 * h_y);
            return df_dy;
        } else if (W == WhichDerivative::xx) {
            double df2_dxdx = (F(x + h_x, y) - 2 * F(x, y) + F(x - h_x, y)) / (h_x * h_x);
            return df2_dxdx;
        } else if (W == WhichDerivative::yy) {
            double df2_dydy = (F(x, y + h_y) - 2 * F(x, y) + F(x, y - h_y)) / (h_y * h_y);
            return df2_dydy;
        } else if (W == WhichDerivative::xy) {
            double df2_dxdy = (F(x + h_x, y + h_y) - F(x + h_x, y - h_y) - F(x - h_x, y + h_y) + F(x - h_x, y - h_y)) /
                              (4 * h_x * h_y);
            return df2_dxdy;
        }
    } else if constexpr (M == DiffMethod::stencil5) {
        if (W == WhichDerivative::x) {
            double df_dx =
                    (F(x - 2 * h_x, y) - F(x + 2 * h_x, y) + 8 * F(x + h_x, y) - 8 * F(x - h_x, y)) / (12 * h_x);
            return df_dx;
        } else if (W == WhichDerivative::y) {
            double df_dy =
                    (F(x, y - 2 * h_y) - F(x, y + 2 * h_y) + 8 * F(x, y + h_y) - 8 * F(x, y - h_y)) / (12 * h_y);
            return df_dy;
        } else if (W == WhichDerivative::xx) {
            double df2_dxdx =
                    (-F(x + 2 * h_x, y) + 16 * F(x + h_x, y) - 30 * F(x, y) + 16 * F(x - h_x, y) - F(x - 2 * h_x, y)) /
                    (12 * h_x * h_x);
            return df2_dxdx;
        } else if (W == WhichDerivative::yy) {
            double df2_dydy =
                    (-F(x, y + 2 * h_y) + 16 * F(x, y + h_y) - 30 * F(x, y) + 16 * F(x, y - h_y) - F(x, y - 2 * h_y)) /
                    (12 * h_y * h_y);
            return df2_dydy;
        } else if (W == WhichDerivative::xy) {
            double df2_dxdy = (F(x + h_x, y + h_y) - F(x + h_x, y - h_y) - F(x - h_x, y + h_y) + F(x - h_x, y - h_y)) /
                              (4 * h_x * h_y);
            return df2_dxdy;
        }
    } else if constexpr (M == DiffMethod::stencil3Extra) {
        if (W == WhichDerivative::x) {
            double df_dx_extra = (n * n * (F(x + h_x / n, y) - F(x - h_x / n, y)) / (2 * h_x / n) -
                                  (F(x + h_x, y) - F(x - h_x, y)) / (2 * h_x)) / (n * n - 1);
            return df_dx_extra;
        } else if (W == WhichDerivative::y) {
            double df_dx_extra = (n * n * (F(x, y + h_y / n) - F(x, y - h_y / n)) / (2 * h_y / n) -
                                  (F(x, y + h_y) - F(x, y - h_y)) / (2 * h_y)) / (n * n - 1);
            return df_dx_extra;
        } else if (W == WhichDerivative::xx) {
            double df2_dxdx_extra = (n * n * (F(x + h_x / n, y) - 2 * F(x, y) + F(x - h_x / n, y)) /
                                     ((h_x * h_x) / (n * n) -
                                      ((F(x + h_x, y) - 2 * F(x, y) + F(x - h_x, y)) / (h_x * h_x))) / (n * n - 1));
            return df2_dxdx_extra;
        } else if (W == WhichDerivative::yy) {
            double df2_dxdx_extra =
                    (n * n * (F(x, y + h_y / n) - 2 * F(x, y) + F(x, y - h_y / n)) / ((h_y * h_y) / (n * n)) -
                     ((F(x, y + h_y) - 2 * F(x, y) + F(x, y - h_y)) / (h_y * h_y))) / (n * n - 1);
            return df2_dxdx_extra;
        } else if (W == WhichDerivative::xy) {
            double df2_dxdy =
                    (n * n * ((F(x + h_x / n, y + h_y / n) - F(x + h_x / n, y - h_y / n) - F(x - h_x / n, y + h_y / n) +
                               F(x - h_x / n, y - h_y / n)) /
                              (4 * h_x / n * h_y / n)) -
                     ((F(x + h_x, y + h_y) - F(x + h_x, y - h_y) - F(x - h_x, y + h_y) + F(x - h_x, y - h_y)) /
                      (4 * h_x * h_y))) / (n * n - 1);

            return df2_dxdy;
        }
    } else if constexpr (M == DiffMethod::stencil5Extra) {
        if (W == WhichDerivative::x) {
            double df_dx =
                    (F(x - 2 * h_x / n, y) - F(x + 2 * h_x / n, y) + 8 * F(x + h_x / n, y) - 8 * F(x - h_x / n, y)) /
                    (12 * h_x / n);
            return df_dx;
        } else if (W == WhichDerivative::y) {
            double df_dy =
                    (F(x, y - 2 * h_y / n) - F(x, y + 2 * h_y / n) + 8 * F(x, y + h_y / n) - 8 * F(x, y - h_y / n)) /
                    (12 * h_y / n);
            return df_dy;
        } else if (W == WhichDerivative::xx) {
            double df2_dxdx =
                    (-F(x + 2 * h_x / n, y) + 16 * F(x + h_x / n, y) - 30 * F(x, y) + 16 * F(x - h_x / n, y) -
                     F(x - 2 * h_x / n, y)) /
                    (12 * h_x / n * h_x / n);
            return df2_dxdx;
        } else if (W == WhichDerivative::yy) {
            double df2_dydy =
                    (-F(x, y + 2 * h_y / n) + 16 * F(x, y + h_y / n) - 30 * F(x, y) + 16 * F(x, y - h_y / n) -
                     F(x, y - 2 * h_y / n)) /
                    (12 * h_y / n * h_y / n);
            return df2_dydy;
        } else if (W == WhichDerivative::xy) {
            double df2_dxdy = (F(x + h_x / n, y + h_y / n) - F(x + h_x / n, y - h_y / n) - F(x - h_x / n, y + h_y / n) +
                               F(x - h_x / n, y - h_y / n)) /
                              (4 * h_x / n * h_y / n);
            return df2_dxdy;
        }
    } else if constexpr (M == DiffMethod::FwdAAD) {
        // TODO
    }
}
