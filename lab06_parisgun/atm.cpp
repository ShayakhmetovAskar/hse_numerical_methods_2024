#include "consts.hpp"
#include "atm.hpp"
#include <cmath>
#include <iostream>
#include <vector>

namespace ADAAI {
    std::pair<double, double> SAM::operator()(double h) {
        double Th;
        double pressure;
        double density;
        double A;
        if (h >= h0 && h <= h1) {
            Th = T0 - r0 * (h - h0);
            pressure = p0 * std::exp(g * std::log(1 - (r0 * (h - h0)) / T0) / (R_air * r0));
            density = pressure / (R_air * (T0 - r0 * (h - h0)));
            A = std::sqrt(pressure / density);
        }
        else if (h > h1 && h <= h2) {
            Th = T1 - r1 * (h - h1);
            pressure = p1 * std::exp(-g * (h - h1) / (R_air * T1));
            density = pressure / (R_air * (T1 - r1 * (h - h1)));
            A = std::sqrt(pressure / density);
        }
        else if (h > h2 && h <= h3) {
            Th = T2 - r2 * (h - h2);
            pressure = p2 * std::exp(g * std::log(1 - (r2 * (h - h2)) / T2) / (R_air * r2));
            density = pressure / (R_air * (T2 - r2 * (h - h2)));
            A = std::sqrt(pressure / density);
        }
        else {
            Th = T3 - r3 * (h - h3);
            pressure = p3 * std::exp(g * std::log(1 - (r3 * (h - h3)) / T3) / (R_air * r3));
            density = pressure / (R_air * (T3 - r3 * (h - h3)));
            A = std::sqrt(pressure / density);
        }
        return std::make_pair(pressure, density);
    }

    double SAM::get_pressure(double h) {
        double Th;
        double pressure;
        if (h >= h0 && h <= h1) {
            Th = T0 - r0 * (h - h0);
            pressure = p0 * std::exp(g * std::log(1 - (r0 * (h - h0)) / T0) / (R_air * r0));
        }
        else if (h > h1 && h <= h2) {
            Th = T1 - r1 * (h - h1);
            pressure = p1 * std::exp(-g * (h - h1) / (R_air * T1));
        }
        else if (h > h2 && h <= h3) {
            Th = T2 - r2 * (h - h2);
            pressure = p2 * std::exp(g * std::log(1 - (r2 * (h - h2)) / T2) / (R_air * r2));
        }
        else {
            Th = T3 - r3 * (h - h3);
            pressure = p3 * std::exp(g * std::log(1 - (r3 * (h - h3)) / T3) / (R_air * r3));
        }
        return pressure;
    }

    double SAM::get_density(double h) {
        double pressure = get_pressure(h);
        double density;
        if (h >= h0 && h <= h1) {
            density = pressure / (R_air * (T0 - r0 * (h - h0)));
        }
        else if (h > h1 && h <= h2) {
            density = pressure / (R_air * (T1 - r1 * (h - h1)));
        }
        else if (h > h2 && h <= h3) {
            density = pressure / (R_air * (T2 - r2 * (h - h2)));
        }
        else {
            density = pressure / (R_air * (T3 - r3 * (h - h3)));
        }
        return density;
    }

    double SAM::get_A(double h) {
        double pressure = get_pressure(h);
        double density = get_density(h);
        double A = std::sqrt(pressure / density);
        return A;
    }
}

int main() {
    ADAAI::SAM atmosphereModel;

    double altitude;
    std::cout << "Введите высоту (в метрах): ";
    std::cin >> altitude;

    std::pair<double, double> result = atmosphereModel(altitude);
    std::cout << "Давление на высоте " << altitude << " м: " << result.first << " Па" << std::endl;
    std::cout << "Плотность на высоте " << altitude << " м: " << result.second << " кг/м^3" << std::endl;
    std::cout << '\n';
    std::cout << "Давление на высоте " << altitude << " м: " << ADAAI::SAM::get_pressure(altitude) << " Па" << std::endl;
    std::cout << "Плотность на высоте " << altitude << " м: " << ADAAI::SAM::get_density(altitude) << " Па" << std::endl;
    std::cout << "Скорость звука на высоте " << altitude << " м: " << ADAAI::SAM::get_A(altitude) << " Па" << std::endl;

    return 0;
}
