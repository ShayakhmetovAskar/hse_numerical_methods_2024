#include "differentiation.hpp"
#include <cmath>

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

    class AAD22 {

    public: // TODO accessors
        double m_val; // value of func or arg
        double m_d1[2]; //df/dx, df/dy in some point (x,y) where value was computed
        // to remember derivative
        double m_d2[3]; //d^2f/dx^2, d^2f/dy^2, d^f/dxdy in the same point
    public:
        AAD22() = delete; // makes no sense
        // construct AAD22 from a const
        constexpr AAD22(double c) :
                m_val(c),
                m_d1{0, 0},
                m_d2{0, 0, 0} {};
    public:
        constexpr AAD22(int i, double v) : // i = 0 for x, i = 1 for y
                m_val(v),
                m_d1{(i == 0) ? 1.0 : 0.0, (i == 0 ? 0.0 : 1.0)},
                m_d2{0, 0, 0} {};
    public:
        constexpr static AAD22 X(double v) {
            return {0, v};
        }

        constexpr static AAD22 Y(double v) {
            return {1, v};
        }

        // unary operator overloading
        AAD22 &operator+() {
            return *this;
        }

        AAD22 operator-() {
            m_val = -m_val;
            m_d1[0] = -m_d1[0];
            m_d1[1] = -m_d1[1];
            m_d2[0] = -m_d2[0];
            m_d2[1] = -m_d2[1];
            m_d2[2] = -m_d2[2];
            return *this;
        }

        // binary operator overloading
        AAD22 operator+(const AAD22 &right) {
            AAD22 res = *this;
            res.m_val += right.m_val;
            res.m_d1[0] += right.m_d1[0];
            res.m_d1[1] += right.m_d1[1];
            res.m_d2[0] += right.m_d2[0];
            res.m_d2[1] += right.m_d2[1];
            res.m_d2[2] += right.m_d2[2];
            return res;
        }

        AAD22 operator-(const AAD22 &right) {
            AAD22 res = *this;
            res.m_val -= right.m_val;
            res.m_d1[0] -= right.m_d1[0];
            res.m_d1[1] -= right.m_d1[1];
            res.m_d2[0] -= right.m_d2[0];
            res.m_d2[1] -= right.m_d2[1];
            res.m_d2[2] -= right.m_d2[2];
            return res;
        }

        AAD22 operator*(const AAD22 &right) {
            AAD22 res = *this;
            res.m_val *= right.m_val;
            res.m_d1[0] = right.m_val * this->m_d1[0] + this->m_val * right.m_d1[0]; //(uv)'=u'v+v'u
            res.m_d1[1] = right.m_val * this->m_d1[1] + this->m_val * right.m_d1[1];
            res.m_d2[0] = this->m_d2[0] * right.m_val + 2 * this->m_d1[0] * right.m_d1[0] + this->m_val * right.m_d2[0];
            res.m_d2[1] = this->m_d2[1] * right.m_val + 2 * this->m_d1[1] * right.m_d1[1] + this->m_val * right.m_d2[1];
            res.m_d2[2] = this->m_d2[2] * right.m_val + this->m_d1[0] * right.m_d1[1] + this->m_val * right.m_d2[2];
            // (fg)''=f''g+2f'g'+fg''
            // TODO needs a check
            return res;
            // [v1, d1, d2] * [v2, g1, g2]
            // [v1*v2, d1v2 + g1v1, d2v2 + g2v1]
        }

        AAD22 operator/(const AAD22 &right) {
            AAD22 res = *this;
            res.m_val /= right.m_val;
            res.m_d1[0] = (this->m_d1[0] * right.m_val - this->m_val * right.m_d1[0]) / (right.m_val * right.m_val);
            res.m_d1[1] = (this->m_d1[1] * right.m_val - this->m_val * right.m_d1[1]) / (right.m_val * right.m_val);
            res.m_d2[0] =
                    (this->m_d2[0] * right.m_val - 2 * this->m_d1[0] * right.m_d1[0] - this->m_val * right.m_d2[0]) /
                    (right.m_val *
                     right.m_val);
            res.m_d2[1] =
                    (this->m_d2[1] * right.m_val - 2 * this->m_d1[1] * right.m_d1[1] - this->m_val * right.m_d2[1]) /
                    (right.m_val *
                     right.m_val);;
            res.m_d2[2] =
                    (this->m_d2[2] * right.m_val - this->m_d1[0] * right.m_d1[1] - this->m_val * right.m_d2[2]) /
                    (right.m_val *
                     right.m_val);;
            // TODO needs a check
            return res;
        }

        AAD22 operator+=(AAD22 const &right) {
            AAD22 res = *this;
            res = res + right;
            return res;
        }

        AAD22 operator-=(AAD22 const &right) {
            AAD22 res = *this;
            res = res - right;
            return res;
        }

        AAD22 operator*=(AAD22 const &right) {
            AAD22 res = *this;
            res = res * right;
            return res;
        }

        AAD22 operator/=(AAD22 const &right) {
            AAD22 res = *this;
            res = res / right;
            return res;
        }

        AAD22 sin() const {
            AAD22 elem = *this;
            elem.m_val = std::sin(this->m_val);
            elem.m_d1[0] = std::cos(this->m_val) * this->m_d1[0];
            elem.m_d1[1] = std::cos(this->m_val) * this->m_d1[1];
            elem.m_d2[0] =
                    std::cos(this->m_val) * this->m_d2[0] - std::sin(this->m_val) * this->m_d1[0] * this->m_d1[0];
            elem.m_d2[1] =
                    std::cos(this->m_val) * this->m_d2[1] - std::sin(this->m_val) * this->m_d1[1] * this->m_d1[1];
            elem.m_d2[2] =
                    -std::sin(this->m_val) * this->m_d1[0] * this->m_d1[1] + std::cos(this->m_val) * this->m_d2[2];
            return elem;
        }

        AAD22 cos() const {
            AAD22 elem = *this;
            elem.m_val = std::cos(this->m_val);
            elem.m_d1[0] = -std::sin(this->m_val) * this->m_d1[0];
            elem.m_d1[1] = -std::sin(this->m_val) * this->m_d1[1];
            elem.m_d2[0] = -std::cos(this->m_val) * this->m_d1[0] * this->m_d1[0] -
                           this->m_d2[0] * this->m_d2[0] * std::sin(this->m_val);
            elem.m_d2[1] = -std::cos(this->m_val) * this->m_d1[1] * this->m_d1[1] -
                           this->m_d2[1] * this->m_d2[1] * std::sin(this->m_val);
            elem.m_d2[2] =
                    -std::cos(this->m_val) * this->m_d1[0] * this->m_d1[1] - this->m_d2[2] * std::sin(this->m_val);
            return elem;
        }

        AAD22 exp() const {
            AAD22 elem = *this;
            elem.m_val = std::exp(this->m_val);
            elem.m_d1[0] = std::exp(this->m_val) * this->m_d1[0];
            elem.m_d1[1] = std::exp(this->m_val) * this->m_d1[1];
            elem.m_d2[0] =
                    std::exp(this->m_val) * this->m_d1[0] * this->m_d1[0] + std::exp(this->m_val) * this->m_d2[0];
            elem.m_d2[1] =
                    std::exp(this->m_val) * this->m_d1[1] * this->m_d1[1] + std::exp(this->m_val) * this->m_d2[1];
            elem.m_d2[2] =
                    std::exp(this->m_val) * this->m_d1[0] * this->m_d1[1] + std::exp(this->m_val) * this->m_d2[2];
            return elem;
        }
    };

    AAD22 sin(AAD22 const &v) {
        return v.sin();
    }

    AAD22 cos(AAD22 const &v) {
        return v.cos();
    }

    AAD22 exp(AAD22 const &v) {
        return v.exp();
    }

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
                double df2_dxdy =
                        (F(x + h_x, y + h_y) - F(x + h_x, y - h_y) - F(x - h_x, y + h_y) + F(x - h_x, y - h_y)) /
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
                        (-F(x + 2 * h_x, y) + 16 * F(x + h_x, y) - 30 * F(x, y) + 16 * F(x - h_x, y) -
                         F(x - 2 * h_x, y)) /
                        (12 * h_x * h_x);
                return df2_dxdx;
            } else if (W == WhichDerivative::yy) {
                double df2_dydy =
                        (-F(x, y + 2 * h_y) + 16 * F(x, y + h_y) - 30 * F(x, y) + 16 * F(x, y - h_y) -
                         F(x, y - 2 * h_y)) /
                        (12 * h_y * h_y);
                return df2_dydy;
            } else if (W == WhichDerivative::xy) {
                double df2_dxdy =
                        (F(x + h_x, y + h_y) - F(x + h_x, y - h_y) - F(x - h_x, y + h_y) + F(x - h_x, y - h_y)) /
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
                                         ((h_x * h_x) / (n * n)) -
                                         ((F(x + h_x, y) - 2 * F(x, y) + F(x - h_x, y)) / (h_x * h_x))) / (n * n - 1);
                return df2_dxdx_extra;
            } else if (W == WhichDerivative::yy) {
                double df2_dxdx_extra =
                        (n * n * (F(x, y + h_y / n) - 2 * F(x, y) + F(x, y - h_y / n)) / ((h_y * h_y) / (n * n)) -
                         ((F(x, y + h_y) - 2 * F(x, y) + F(x, y - h_y)) / (h_y * h_y))) / (n * n - 1);
                return df2_dxdx_extra;
            } else if (W == WhichDerivative::xy) {
                double df2_dxdy =
                        (n * n *
                         ((F(x + h_x / n, y + h_y / n) - F(x + h_x / n, y - h_y / n) - F(x - h_x / n, y + h_y / n) +
                           F(x - h_x / n, y - h_y / n)) /
                          (4 * h_x / n * h_y / n)) -
                         ((F(x + h_x, y + h_y) - F(x + h_x, y - h_y) - F(x - h_x, y + h_y) + F(x - h_x, y - h_y)) /
                          (4 * h_x * h_y))) / (n * n - 1);

                return df2_dxdy;
            }
        } else if constexpr (M == DiffMethod::stencil5Extra) {
            if (W == WhichDerivative::x) {
                double df_dx =
                        (F(x - 2 * h_x / n, y) - F(x + 2 * h_x / n, y) + 8 * F(x + h_x / n, y) -
                         8 * F(x - h_x / n, y)) /
                        (12 * h_x / n);
                return df_dx;
            } else if (W == WhichDerivative::y) {
                double df_dy =
                        (F(x, y - 2 * h_y / n) - F(x, y + 2 * h_y / n) + 8 * F(x, y + h_y / n) -
                         8 * F(x, y - h_y / n)) /
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
                double df2_dxdy =
                        (F(x + h_x / n, y + h_y / n) - F(x + h_x / n, y - h_y / n) - F(x - h_x / n, y + h_y / n) +
                         F(x - h_x / n, y - h_y / n)) /
                        (4 * h_x / n * h_y / n);
                return df2_dxdy;
            }
        } else if constexpr (M == DiffMethod::FwdAAD) {
            AAD22 X = AAD22::X(x);
            AAD22 Y = AAD22::Y(y);
            AAD22 Res = F(X, Y); // F - our callable
            // Call the corresponding accessor depending on the WhichDerivative template arg
            // TODO
            if (W == WhichDerivative::x) {
                return Res.m_d1[0];
            } else if (W == WhichDerivative::y) {
                return Res.m_d1[1];
            } else if (W == WhichDerivative::xx) {
                return Res.m_d2[0];
            } else if (W == WhichDerivative::yy) {
                return Res.m_d2[1];
            } else if (W == WhichDerivative::xy) {
                return Res.m_d2[2];
            }
        }
    }
}
