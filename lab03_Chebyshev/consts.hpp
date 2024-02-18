#ifndef ADAAI_CONSTS_HPP
#define ADAAI_CONSTS_HPP

#include <cmath>
#include <climits>

namespace ADAAI {
    template<typename F>
    constexpr inline F c_Eps = std::numeric_limits<F>::epsilon();
}
#endif