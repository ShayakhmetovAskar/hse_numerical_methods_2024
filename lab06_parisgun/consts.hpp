#ifndef ADAAI_CONSTS_HPP
#define ADAAI_CONSTS_HPP

#include <cmath>
#include <climits>

namespace ADAAI {
    constexpr double h0 = 0.0;
    constexpr double h1 = 11000.0;
    constexpr double h2 = 20000.0;
    constexpr double h3 = 32000.0;
    constexpr double h4 = 47000.0;

    constexpr double r0 = 0.0065;
    constexpr double r1 = 0.0;
    constexpr double r2 = -0.001;
    constexpr double r3 = -0.0028;

    constexpr double T0 = 288.15; // 0 слой
    constexpr double T1 = 213.65; // 1 слой
    constexpr double T2 = 213.65; // 2 слой
    constexpr double T3 = 225.65; // 3 слой
    constexpr double T4 = 267.65; // 4 слой

    constexpr double g = 9.80655;

    constexpr double R_air = 287.0528;

    constexpr double p0 = 101325.0; // давление в 0 слое
    constexpr double p1 = 22632.38; // давление в 1 слое
    constexpr double p2 = 5367.01; // давление во 2 слое
    constexpr double p3 = 868.02; // давление в 3 слое

    constexpr double mass = 106.0;
    constexpr double diameter = 0.216;
    constexpr double v_0 = 1640.0;
    constexpr double alpha = 55.0;

    constexpr double pi = M_PI;

    double S = pi * diameter * diameter / 4.0;
    double alpha_rad = pi * alpha / 180.0;

    constexpr double omega = 7.2921e-5;
    constexpr double R = 6.371e6;
}
#endif