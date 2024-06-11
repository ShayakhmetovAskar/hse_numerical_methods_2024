#ifndef ADAAI_DIFFERENTIATION_HPP
#define ADAAI_DIFFERENTIATION_HPP

#include <climits>
#include <cmath>
#include <stdexcept>

namespace ADAAI {
enum class DiffMethod : int {
  stencil3,
  stencil3Extra,
  stencil5,
  stencil5Extra,
  FwdAAD
};

enum class WhichDerivative : int {
  x,  // df/dx
  y,  // df/dy
  xx, // d^2f/(dx)^2
  yy, // d^2f/(dy)^2
  xy  // d^2f/dxdy
};
class AAD22 {
private:
  double m_val;   // value of func or arg
  double m_d1[2]; // df/dx, df/dy in some point (x,y) where value was computed
  // to remember derivative
  double m_d2[3]; // d^2f/dx^2, d^2f/dy^2, d^f/dxdy in the same point

public:
  double get_val() const;
  void set_val(double val);
  const double *get_d1() const;
  void set_d1(double d1x, double d1y);
  const double *get_d2() const;
  void set_d2(double d2xx, double d2yy, double d2xy);

  constexpr explicit AAD22(double c);
  constexpr AAD22(int i, double v);

  constexpr static AAD22 X(double v);
  constexpr static AAD22 Y(double v);

  AAD22 &operator+();
  AAD22 operator-();

  AAD22 operator+(const AAD22 &right);
  AAD22 operator-(const AAD22 &right);
  AAD22 operator*(const AAD22 &right);
  AAD22 operator/(const AAD22 &right);

  AAD22 operator+=(AAD22 const &right);
  AAD22 operator-=(AAD22 const &right);
  AAD22 operator*=(AAD22 const &right);
  AAD22 operator/=(AAD22 const &right);

  AAD22 sin() const;
  AAD22 cos() const;
  AAD22 exp() const;
};

AAD22 sin(AAD22 const &v);
AAD22 cos(AAD22 const &v);
AAD22 exp(AAD22 const &v);

template <WhichDerivative W, DiffMethod M, typename Callable>
double Differentiator(Callable const &F, double x, double y);

} // namespace ADAAI
#endif // ADAAI_DIFFERENTIATION_HPP