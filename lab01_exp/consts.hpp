#pragma once

#include <type_traits>
#include <cmath>
#include <cfloat>

namespace ADAAI {
    template<typename F>
    constexpr inline F ln2;

    template<>
    constexpr inline float ln2<float> = (float) M_LOG2E;
    template<>
    constexpr inline double ln2<double> = M_LOG2E;
    template<>
    constexpr inline long double ln2<long double> = (long double)M_LOG2E;

    template<typename F>
    constexpr inline F eps;
    template<>
    constexpr inline float eps<float> = FLT_EPSILON;
    template<>
    constexpr inline double eps<double> = DBL_EPSILON;
    template<>
    constexpr inline long double eps<long double> = LDBL_EPSILON;

    template<typename F>
    constexpr inline F sqrt2;
    template<>
    constexpr inline float sqrt2<float> = (float) M_SQRT2;
    template<>
    constexpr inline float sqrt2<double> = (double) M_SQRT2;
    template<>
    constexpr inline float sqrt2<long double> = (long double) M_SQRT2;

    template<typename F>
    constexpr inline F inf;
    template<>
    constexpr inline float inf<float> = HUGE_VALF;
    template<>
    constexpr inline float inf<double> = HUGE_VAL;
    template<>
    constexpr inline float inf<long double> = HUGE_VALL;
}