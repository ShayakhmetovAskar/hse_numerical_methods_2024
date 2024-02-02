#ifndef ADAAI_CONSTS_HPP
#define ADAAI_CONSTS_HPP

#include <cmath>
#include <climits>

namespace ADAAI {
    template<typename F>
    constexpr inline F ln2;

    template<>
    constexpr inline float ln2<float> = 1 / M_LOG2Ef;
    template<>
    constexpr inline double ln2<double> = 1 / M_LOG2E;
    template<>
    constexpr inline long double ln2<long double> = 1 / M_LOG2El;

    template<typename F>
    constexpr inline F eps = std::numeric_limits<F>::epsilon();

    template<typename F>
    constexpr inline F sqrt2;
    template<>
    constexpr inline float sqrt2<float> = M_SQRT2f;
    template<>
    constexpr inline double sqrt2<double> = M_SQRT2;
    template<>
    constexpr inline long double sqrt2<long double> = M_SQRT2l;

    template<typename F>
    constexpr inline F inf;
    template<>
    constexpr inline float inf<float> = HUGE_VALF;
    template<>
    constexpr inline float inf<double> = HUGE_VAL;
    template<>
    constexpr inline float inf<long double> = HUGE_VALL;
}
#endif