#ifndef ADAAI_CONSTS_HPP
#define ADAAI_CONSTS_HPP

#include <climits>
#include <cmath>

namespace ADAAI {
template <typename F> constexpr inline F c_Ln2;

template <> constexpr inline float c_Ln2<float> = 1 / (float)M_LOG2E;
template <> constexpr inline double c_Ln2<double> = 1 / M_LOG2E;
template <>
constexpr inline long double c_Ln2<long double> = 1 / (long double)M_LOG2E;

template <typename F>
constexpr inline F c_Eps = std::numeric_limits<F>::epsilon();

template <typename F> constexpr inline F c_Inf;
template <> constexpr inline float c_Inf<float> = HUGE_VALF;
template <> constexpr inline float c_Inf<double> = HUGE_VAL;
template <> constexpr inline float c_Inf<long double> = HUGE_VALL;
} // namespace ADAAI
#endif