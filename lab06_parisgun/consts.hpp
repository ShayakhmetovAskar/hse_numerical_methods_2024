#ifndef ADAAI_CONSTS_HPP
#define ADAAI_CONSTS_HPP

#include <cmath>
#include <climits>

namespace ADAAI {
    const double h0 = 0.0;
    const double h1 = 11000.0;
    const double h2 = 20000.0;
    const double h3 = 32000.0;
    const double h4 = 47000.0;

    const double r0 = 0.0065;
    const double r1 = 0.0;
    const double r2 = -0.001;
    const double r3 = -0.0028;

    const double T0 = 288.15; // 0 слой
    const double T1 = 213.65; // 1 слой
    const double T2 = 213.65; // 2 слой
    const double T3 = 225.65; // 3 слой
    const double T4 = 267.65; // 4 слой

    const double g = 9.80655;

    const double R_air = 287.0528;

    const double p0 = 101325.0; // давление в 0 слое
    const double p1 = 22632.38; // давление в 1 слое
    const double p2 = 5367.01; // давление во 2 слое
    const double p3 = 868.02; // давление в 3 слое

    const double mass = 106.0;
    const double diameter = 0.216;
    const double v_0 = 1640.0;
    const double alpha = 55.0;

    const double pi = M_PI;

    double S = pi * diameter * diameter / 4.0;
    double alpha_rad = pi * alpha / 180.0;


}
#endif